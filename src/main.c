//
//	Author: Alexander B. Ringheim
//	Project: FEM in Plasticity
//	Date of creation: 16.10.2024
//

// System includes
#include <stdio.h>

// Local includes
#include "croutReduction.h"
#include "elements.h"
#include "forceVector.h"
#include "jacobian.h"
#include "returnMapping.h"
#include "stiffnessMatrix.h"
#include "skyline.h"
#include "solveCroutSkyline.h"
#include "shapeFunctions.h"
#include "vonMisesYield.h"

int main (char *args)
{
    // The entry point of the program.

    /*    
    Step 1: Import and initialize Data
        - Load input data: Node coords, Element connectivity, Boundary conditions, Applied load vectors, Material properties (E, v, sigma_yield_0, H)
        - Compute and store shape functions and its derivatives in the reference space (xi, eta) (static values)
        - Define Gauss Point locations and weights for numerical integration (static values)

    Step 2: Initialize Gauss Point internal variables
        - For each element and Gauss Point:
            - epsilon_p = 0
            - epsilon_bar_p = 0
            - sigma = 0
        - These variables will be updated as plasticity accumulates per load step
    
    From here, everything is within the iterator, k
    Step 3: Newton-Raphson loop starts
        - Apply load increment (f_ext + delta f_ext)
        - Get the initial displacement "guess" (not actually a guess. First is the initial value, often u0 = 0, and from each step it will be incremented from the solution of the system).

    Step 4: Loop over elements & Gauss points
        - Extract the nodal displacements for the element from u, -> u_e
        - Compute the Jacobian and B-matrix
        - Compute the strain as epsilon = B * u_e
        - Compute the trial stress as sigma_trial = D_e (epsilon - epsilon_p)
        - From the trial stress, calculate the trial deviatoric stress and the von Mises equivalent stress
        - Check the trial stress against the yield function f = sigma_eq - sigma_y, where sigma_y = sigma_y0 + H * epsilon_bar_p
        - If f > 0 -> Plasticity (compute delta gamma, n, and D_ep), if f <= 0 -> elasticity (use D_e)
            - if f > 0 -> Compute delta gamma (plastic multiplier) = f_yield / (3 G + H)
            - G = E / (2(1+v))
            - Correct stress on Gauss point as sigma = sigma_trial - 2 G * delta gamma * n
            - Make trial epsilon_p = epsilon_p_stored + delta gamma * n
            - Make trial epsilon_bar_p = epsilon_bar_stored_p + delta gamma * sqrt(2/3)
                - THESE 3 VARIABLES, SIGMA, TRIAL EPSILON_P AND TRIAL_EPSILON_BAR_P ARE ONLY COMMITED AFTER CONVERGENCE CHECK.
    
    Step 5: Compute Local Element Quantities (could be included in step 4 since it is still in each element loop)
        - Compute K_e = sum(gp) Btrans D(e / ep determined in step 4) B * detJ * weight(gp)
        - Compute f_int_e = sum(gp) Btrans * sigma * detJ * weight(gp)
    
    Step 6: Assemble global system
        - Assemble global stiffness matrix, K
        - Assemble global internal force vector

    Step 7: Solve system and update displacement
        - r^(k) = f_ext^(k) - f_int^(k)
        - K^(k) * delta u^(k) = r^(k)
        - Increment u^(k+1) = u^(k) + delta u^(k)

    Step 8: Convergence check for the current step
        - Check that abs(r^(k)) < tolerance or abs(delta u^(k)) is small
        - if converged ->
            - Commit values for epsilon_p, epsilon_bar_p and sigma for each Gauss Point
            - Add a load increment and go to the next Newton-Raphson step. (Go to step 3 and add to f_ext).
        - else 
            - Go back to step 3 but do not add a load increment if not converged!
    
    Step 9: Once the entire load is applied, and the system has converged
        - Run post-processing to extract results in useful formats.
    */
    printf("Running Finite Elements Simulation on Plasticity Equations! \n");



    // DEFINING VARIABLES AND ALLOCATING MEMORY FOR THE PROCESSES



    // Material properties
    double E; // Young's Modulus
    double v; // Poisson's ratio
    double sigmaYieldInitial; // Initial yield stress
    double H; // Hardening modulus
    double G = shearModulus(E, v);
    
    // Constant values for the analysis
    int gp = 4; // Gauss Points per element
    int nelements = 100; // Number of elements in the mesh.
    int nnodesElement = 4; // Nodes per element
    int nnodes = 200; // Number of nodes in the mesh.
    int DOF = 2; // Degrees of freedom per node
    int dim = 2; // Dimensions of the system
    int Km = nnodes * DOF; // Number of rows / columns of the global stiffness matrix.
    double *K;
    skylineMatrix *Kskyline;

    int maxNodesBeforeSkyline = 10000;
    int skylineSolver;

    if (Km * Km > maxNodesBeforeSkyline)
    {
        // Skyline
        //Kskyline = initSkylineMatrix(Km);
        skylineSolver = 1;
    }
    else
    {
        // Normal matrix
        skylineSolver = 0;
        K = malloc(Km * Km * sizeof(double));
        initGlobalStiffnessMatrix(K, Km); // Set all values in the matrix to 0.
    }

    // Allocate memory for the global internal force vector, f_int, the global external force vector, f_ext, and the residual, r.
    double *Fint = malloc(Km * sizeof(double));
    double *Fext = malloc(Km * sizeof(double));
    double *r = malloc(Km * sizeof(double));

    double *Fstep = malloc(Km * sizeof(double)); // This is the vector that will hold the incremented load step.
    for (int i = 0; i < Km; i++)
    {
        *(Fstep + i) = 0; // Set the intial values to 0.
    }

    // Determine whether to use normal matrix or skyline matrix.

    // Allocate memory for the 'element' and 'node' structs.
    quadElement *element = malloc(nelements * sizeof(quadElement));

    // Allocate memory for the global-, increment- and element displacement vector.
    double *u = malloc(nnodes * DOF * sizeof(double));
    double *du = malloc(nnodes * DOF * sizeof(double));
    double *ue = malloc(nnodesElement * DOF * sizeof(double));

    // Get the shape functions for 2D Quad elements
    shapeFunctions2D *sf = malloc(gp * sizeof(shapeFunctions2D));

    // Allocate memory for the Jacobian matrix, its determinant and the inverse.
    int j_m = 2; // Representing the size, m, in a m*m Jacobian matrix.
    double *J = malloc(j_m * j_m * sizeof(double)); // This is initialized with values in the Jacobian function.
    double *detJ = malloc(sizeof(double)); // Store the determinant value of the Jacobian here.
    double *Jinv = malloc (j_m * j_m * sizeof(double)); // This is also initialized in the Jacobian function.

    // Allocate memory for the B matrix, B transpose, D and the element stiffness matrix.
    int nlocnode = 4; // Number of local nodes per element.
    int Bm = 3; // Rows of the B matrix.
    int Bn = nlocnode * 2; // Each node adds 2 columns to the B-matrix.
    int Dm = 3;
    int Dn = 3;
    double *B = malloc(Bn * Bm * sizeof(double)); 
    double *Btrans = malloc(Bn * Bm * sizeof(double));
    double *D = malloc(Dm * Dn * sizeof(double));

    // Allocate memory for the element stiffness matrix, Ke
    int Kem = nnodesElement * DOF;
    int Ken = Kem;
    double *Ke = malloc(Kem * Ken * sizeof(double));

    // Allocate memory for the element internal force vector, FintE
    double *FintE = malloc(Kem * sizeof(double));
    

    // Variables used in the Newton-Raphson iteration.
    int stepCounter = 0;
    int loadIncrementSteps = 100; // Divide the increments into this many load increments.
    int maxSteps = 100;     // Maximum allowed steps for each converged step.
    int maxLoadSteps = 100; // Maximum allowed steps for each load step.
    int stepConverged = 1; // Check for seeing if convergence has been reached. Set it initially to 1, so that the first load increment will not be 0.
    int fullLoadApplied = 0; // Check for seeing if the entire load is applied. Set initially to 0, and set to 1 only if all the Fload[i] >= Fext[i].
    int solutionFound = 0; // Check for seeing if the solution is found and converged (1), or if the iteration simply ended because the maxSteps was reached.
    int residualThreshold = 0.1; // Threshold for the absolute value of the residual for the convergence criteria.

    // Stresses and strains
    double *epsilon = malloc(nnodesElement * DOF * sizeof(double));

    double *sigmaTrial = malloc(Dn * sizeof(double)); // Trial stress
    double *sDeviatoric = malloc(Dn * sizeof(double)); // Deviatoric stress
    double *nUnitDeviatoric = malloc(Dn * sizeof(double));
    double *sigma = malloc(Dn * sizeof(double)); // Corrected stress for Gauss Point
    double sigmaEq; // von Mises equivalent stress
    double sigmaYield; // Corrected yield stress
    double f; // von Mises yield function result
    double deltaGamma;
    double *trialEpsilonP = malloc(Dn * sizeof(double));
    double trialEpsilonBarP;


    // END OF DEFINING VARIABLES AND ALLOCATING MEMORY


    // Check how many nodes and degrees of freedom there are. If less than threshold, use normal matrix. If over threshold, use Skyline.
    // Allocate memory for the global stiffness matrix, (using Skyline).


    /*
        This marks the actual start of the steps of the Finite Element code...
    */
    // Initialize -> Material Properties -> Stress and Strain Variables -> Initial Conditions and Read Input Data (Geometry, Boundary Conditions & Loading)
    

    
    // Gets the shape functions for a 2D quad element. Needs only be ran once as long as all the elements are the same.
    quadShapeFunctions(sf); 
    


    // Newton-Raphson Iterator
    printf("Starting Newton-Raphson iteration.");
    for (int k = 0; k < maxSteps; k++)
    {
        // The following logic is complex, but it does 4 things.
        // 1: Commit the trial values in the Gauss Point to be true values.
        // 2: Check if all the load is applied AND the last step converged. If all values are past the mark, break out of the loop. 
        // 3, 4: If not, apply load increment AND reset the displacement increment vector, du to 0.
        if (stepConverged)
        {
            commitTrialValuesAtGaussPoints(element, nelements); // Commit trial values to commited values.
            fullLoadApplied = 1; // Try this, if any of the values are not done, then this will be set back to 0 sometime in the for-loop.
            for (int i = 0; i < Km; i++)
            {
                if ( *(Fstep + i) >= *(Fext + i) )
                {
                    fullLoadApplied = 0;
                }
                else
                {
                    *(Fstep + i) += *(Fext + i) / loadIncrementSteps;
                    *(du + i) = 0;
                }
            }
            if (fullLoadApplied)
            {
                printf("Full loads are applied and solution converged on step %d. Proceeding to post-processing . . . \n\n", stepCounter);
                solutionFound = 1; // The solution has been found.
                break; // Finish the loop.
            }
        }

        printf("Running step: %d . . . \n", k);

        // Init / reset global stiffness matrix and global internal force vetor for the current step
        if (skylineSolver)
        {
            // Init / reset skyline matrix
            resetSkylineMatrix(Kskyline);
        }
        else
        {
            // Init / reset regular matrix
            initGlobalStiffnessMatrix(K, Km);
        }
        initGlobalInternalForceVector(Fint, Km);

        // Set convergence check to 0 (false) before the load steps start.
        stepConverged = 0;

        // Start load step
        for (int l = 0; l < maxLoadSteps; l++)
        {
            printf("\tStep: %d, Load Step: %d . . . Running . . . \n\n", k, l);
            stepCounter += 1;

            // Loop through all elements
            for (int e = 0; e < nelements; e++)
            {
                printf("\t\tStep: %d, Load Step: %d . . . Analysing element #%d . . . \n", k, l, e);

                // Init / reset K_e and F_int_e for the current element
                initElementStiffnessMatrix(Ke, Kem);
                initElementInternalForceVector(FintE, Kem);

                // Get the displacements for the element's nodes, u_e, and store them. One value for each DOF.
                for (int i = 0; i < nnodesElement; i++)
                {
                    for (int d = 0; d < DOF; d++) // Loop over each DOF
                    {
                        *(ue + i * DOF + d) = *(u + element[e].nodeids[i] * DOF + d); // Store the displacement from the node with this id.
                    }
                }

                // Begin looping over the Gauss Points of the elements, each loop is for the i-th Gauss Point.
                for (int i = 0; i < element->gp; i++)
                {
                    printf("\t\t\tElement #%d . . . Gauss Point #%d", e, i);

                    // Find the Jacobian, B matrix, and elastic D_e matrix for this Gauss Point.
                    quadJacobian(J, Jinv, j_m, detJ, sf[i].Ni, sf[i].NiPxi, sf[i].NiPeta, element[e].coords, DOF);
                    quadBMatrix(B, Btrans, Jinv, sf[i].NiPxi, sf[i].NiPeta);
                    elasticDMatrixPlaneStrain(D, E, v);

                    // Get the strain for this Gauss Point with the displacement and B matrix.
                    displacementStrain(epsilon, B, ue, nnodesElement, DOF);

                    // Make trial stress
                    trialStress(sigmaTrial, D, epsilon, element[e].epsilonP, Dn, i);

                    // Find deviatoric stress, von Mises equivalent stress, yield stress adjusted for existing plastic strain, and find von Mises Yield Function
                    deviatoricStress2D(sDeviatoric, sigmaTrial);
                    sigmaEq = vonMisesEquivalentStress2D(sDeviatoric);
                    sigmaYield = plasticCorrectedYieldStress(sigmaYieldInitial, H, element[e].epsilonBarP[i]);

                    // Check von Mises Yield Function for yielding of the node. If elastic, use elastic material stiffness matrix. If plastic, run return mapping.
                    f = vonMisesYieldFunction(sigmaEq, sigmaYield);

                    if (f > 0) // f > 0 -> yield
                    {
                        // Return mapping processes

                        // Unit deviatoric stress, n
                        unitDeviatoricStress(nUnitDeviatoric, sDeviatoric, Dn);

                        // Elasto-plastic material stiffness matrix, D_ep
                        elastoPlasticDMatrix(D, nUnitDeviatoric, H, Dn);

                        // Plastic multiplier
                        deltaGamma = plasticMultiplier(f, G, H);

                        // Corrected trial stress
                        plasticStressCorrection(sigma, sigmaTrial, nUnitDeviatoric, G, deltaGamma, Dn);

                        // Append the contribution to the element internal force vector, f^e_int, with the corrected stress. The element internal force vector is summed for each element over all the Gauss Points.
                        elementInternalForceVector(FintE, Btrans, sigma, *detJ, sf[i].weight, Kem, Dn);

                        // Trial plastic strain tensor
                        trialPlasticStrain(trialEpsilonP, element[e].epsilonP, nUnitDeviatoric, deltaGamma, Dn, i);

                        // Trial plastic equivalent strain
                        element[e].trialEpsilonBarP[i] = trialEquivalentPlasticStrain(deltaGamma);

                        // Store trial values in element, but they are not to be commited before the load step converges.
                        for (int j = 0; j < Dn; j++)
                        {
                            element[e].sigma[i * element[e].gp + j] = *(sigma + j);
                            element[e].trialEpsilonP[i * element[e].gp + j] = *(trialEpsilonP + j);
                        }
                    }
                    else
                    {
                        // Append the contribution to the element internal force vector, f^e_int, with the trial stress. The element internal force vector is summed for each element over all the Gauss Points.
                        elementInternalForceVector(FintE, Btrans, sigmaTrial, *detJ, sf[i].weight, Kem, Dn);
                    }

                    // Calculate and append the contribution to the element stiffness matrix, Ke, for the current Gauss Point. The element stiffness matrix is summed for each element over all the Gauss Points.
                    elementStiffnessMatrix(Ke, nnodesElement, DOF, B, Btrans, D, *detJ, sf[i].weight);
                }

                // Assemble into global stiffness matrix.
                if (skylineSolver)
                {
                    // Assemble into skyline matrix
                    globalSkylineStiffnessMatrix(Kskyline,Ke, element[e].nodeids, DOF, nnodesElement, Kem);
                }
                {
                    // Assemble into regular matrix
                    globalStiffnessMatrix(K, Ke, element[e].nodeids, DOF, nnodesElement, Km, Kem);
                }
                
                // Assemble into global internal force vector.
                globalInternalForceVector(Fint, FintE, element[e].nodeids, DOF, nnodesElement, Km, Kem);
            }

            // Compute the residual, r, as f_ext - f_int, and check convergence
            stepConverged = 1; // Make this prediction now, and correct it if ANY of the values in the residual are not below the threshold.
            for (int i = 0; i < Km; i++)
            {
                *(r + i) = *(Fstep + i) - *(Fint + i);
                if ( abs( *(r + i) ) > residualThreshold )
                {
                    stepConverged = 0; // If the current residual value is greater than the threshold, the solution is not converged.
                }
            }

            // If not converged, solve the system.
            if (!stepConverged)
            {
                // Solve system for current applied load.
                if (skylineSolver)
                {
                    // Solve skyline matrix
                    croutSkyline(Kskyline);
                }
                else
                {
                    // Solve regular matrix
                    croutReduction(K, Km, du, r);
                }
            }
            else
            {
                printf("\tStep: %d, Load Step: %d . . . CONVERGED . . . \n\n", k, l);
                break;
            }
        }

        printf("\n\nNewton-Raphson iteration has concluded after %d steps! \n\n", stepCounter);
 
    }


    // Post-processing
    if (solutionFound)
    {

    }
    else
    {
        // No solution was found. Print error message.
    }


    // Freeing some variables from memory now that they are no longer needed.
    free(element); // Element structs
    free(sf); // Shape functions
    free(J); // Jacobian matrix.
    free(detJ); // Determinant of the Jacobian matrix.
    free(Jinv); // Inverse matrix of the Jacobian.
    free(B); // B matrix
    free(Btrans); // B transpose
    free(D);
    free(Ke); // Ke. element stiffness matrix

   
    return 0;
}