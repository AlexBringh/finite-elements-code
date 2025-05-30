//
//	Author: Alexander B. Ringheim
//	Project: FEM in Plasticity
//	Date of creation: 16.10.2024
//

// Standard Module includes
#include <stdio.h>
#include <math.h>

// FEM module includes
#include "croutReduction.h"
#include "elements.h"
#include "forceVector.h"
#include "inputData.h"
#include "jacobian.h"
#include "matrixUtils.h"
#include "postProcessing.h"
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
    double G;
    
    // Constant values for the analysis
    int gp; // Gauss Points per element
    int nelements; // Number of elements in the mesh.
    int nnodesElement; // Nodes per element
    int nnodes; // Number of nodes in the mesh.
    int DOF; // Degrees of freedom per node
    int dim; // Dimensions of the system
    int nmaterials;

    node *nodes = NULL;
    material *materials = NULL;
    quadElement *elements = NULL;

    int *uFixed = NULL;
    double *Fext = NULL;

    // Load data
    if (readCSV("sample/data.csv", &nodes, &nnodes, &nnodesElement, &nelements, &gp, &materials, &nmaterials, &uFixed, &Fext)) {
        return 1;
    }
    DOF = nodes[0].dof;
    dim = nodes[0].dim;
    readElements("sample/data.csv", &elements, &nelements, nodes, nnodes, gp, nnodesElement, DOF);
    E = materials[0].E;
    v = materials[0].v;
    H = materials[0].H;
    G = shearModulus(E, v);
    sigmaYieldInitial = materials[0].yield;
    int Km = nnodes * DOF; // Number of rows / columns of the global stiffness matrix.
    int Kem = nnodesElement * DOF;

    printf("System data: NNodes:%d, NElements:%d, NNodesElement:%d, DOF:%d, Dimensions:%d, NMaterials:%d, GloablDOFs:%d \n", nnodes, nelements, nnodesElement, DOF, dim, nmaterials, Km);

    printf("node\t\tu_x_fixed\tu_y_fixed\tf_x_ext\t\tf_y_ext \n");
    for (int i = 0; i < nnodes; i++)
    {
        printf("%d\t\t%d\t\t%d\t\t%.2f\t\t%.2f \n", i, *(uFixed + i*DOF + 0), *(uFixed + i*DOF + 1), *(Fext + i * DOF + 0), *(Fext + i * DOF + 1));
    }

    
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
        K = calloc(Km * Km, sizeof(double));
        initGlobalStiffnessMatrix(K, Km); // Set all values in the matrix to 0.
    }

    // Allocate memory for the global internal force vector, f_int, the global external force vector, f_ext, and the residual, r.
    double *Fint = calloc(Km, sizeof(double));
    double *r = calloc(Km, sizeof(double));
    double *Fstep = calloc(Km, sizeof(double)); // This is the vector that will hold the incremented load step.

    // Allocate memory for the global-, increment- and element displacement vector, as well as the fixed displacement indices.
    double *u = calloc(Km, sizeof(double));
    double *du = calloc(Km,  sizeof(double));
    double *ue = calloc(Kem, sizeof(double));

    // Get the shape functions for 2D Quad elements
    shapeFunctions2D *sf = malloc(gp * sizeof(shapeFunctions2D));

    // Allocate memory for the Jacobian matrix, its determinant and the inverse.
    int Jn = 2; // Representing the size, m, in a m*m Jacobian matrix.
    double *J = calloc(dim * dim, sizeof(double)); // This is initialized with values in the Jacobian function.
    double *detJ = calloc(1, sizeof(double)); // Store the determinant value of the Jacobian here.
    double *Jinv = calloc(dim * dim, sizeof(double)); // This is also initialized in the Jacobian function.

    // Allocate memory for the B matrix, B transpose, D and the element stiffness matrix.
    int Bm = 3; // Rows of the B matrix.
    int Bn = nnodesElement * 2; // Each node adds 2 columns to the B-matrix.
    int Dn = 3;
    double *B = calloc(Bn * Bm, sizeof(double)); 
    double *Btrans = calloc(Bn * Bm, sizeof(double));
    double *D = calloc(Dn * Dn, sizeof(double));

    // Allocate memory for the element stiffness matrix, Ke
    double *Ke = calloc(Kem * Kem, sizeof(double));

    // Allocate memory for the element internal force vector, FintE
    double *FintE = calloc(Kem, sizeof(double));
    

    // Variables used in the Newton-Raphson iteration.
    int stepCounter = 0;
    int loadIncrementSteps = 100; // Number of load increment steps. Percent = 1/loadIncrementSteps applied load at each step. 100 = 1%, 50 = 2%, 25 = 4%, 20 = 5%, 10 = 10%
    int maxLoadSteps = 50; // Maximum allowed steps for each load step.
    int stepConverged = 1; // Check for seeing if convergence has been reached. Set it initially to 1, so that the first load increment will not be 0.
    int fullLoadApplied = 0; // Check for seeing if the entire load is applied. Set initially to 0, and set to 1 only if all the Fload[i] >= Fext[i].
    int solutionFound = 0; // Check for seeing if the solution is found and converged (1), or if the iteration simply ended because the max number of steps was reached.
    double residualThreshold = 0.1; // Threshold for the absolute value of the residual for the convergence criteria.
    double dispIncrThreshold = 1; // Largest displacement increment before solution is thrown.

    // Stresses and strains
    double *epsilon = calloc(nnodesElement * DOF, sizeof(double));
    double *sigmaTrial = calloc(Dn, sizeof(double)); // Trial stress
    double *sDeviatoric = calloc(Dn, sizeof(double)); // Deviatoric stress
    double *nUnitDeviatoric = calloc(Dn, sizeof(double));
    double *sigma = calloc(Dn, sizeof(double)); // Corrected stress for Gauss Point
    double sigmaEq; // von Mises equivalent stress
    double sigmaYield; // Corrected yield stress
    double f; // von Mises yield function result
    double deltaGamma;
    double *trialEpsilonP = calloc(Dn, sizeof(double));
    double trialEpsilonBarP;


    // END OF DEFINING VARIABLES AND ALLOCATING MEMORY


    // Gets the shape functions for a 2D quad element. Needs only be ran once as long as all the elements are the same.
    quadShapeFunctions(sf); 
    printf("Shape function results: \n");
    for (int i = 0; i < gp; i++)
    {
        printf("N%d,  \t\t N%dPxi, \t\tN%dPeta, \t\tWeight:%.2f \n", (i+1), (i+1), (i+1), sf[i].weight);
        for (int j = 0; j < 4; j++)
        {
            printf("%.4f \t\t %.4f \t\t%.4f \n", sf[i].Ni[j], sf[i].NiPxi[j], sf[i].NiPeta[j]);
        }
        printf("\n");
    }

    printf("Settings: Load Increments:%d, Max load steps:%d, Residual Threshold:%.4f, Displacement Increment Threshold:%.4f \n\n", loadIncrementSteps, maxLoadSteps, residualThreshold, dispIncrThreshold);

    // Newton-Raphson Iterator
    printf("Starting Newton-Raphson iteration. \n\n");
    for (int k = 0; k < loadIncrementSteps + 1; k++)
    {
        // The following logic is complex, but it does 3 things.
        // 1: Check if all the load is applied AND the last step converged. If all values are past the mark, break out of the loop. 
        // 2, 3: If not, apply load increment AND reset the displacement increment vector, du to 0.
        if (stepConverged)
        {
            fullLoadApplied = 1; // Try this, if any of the values are not done, then this will be set back to 0 sometime in the for-loop.
            for (int i = 0; i < Km; i++)
            {
                if ( *(Fstep + i) < *(Fext + i) )
                {
                    fullLoadApplied = 0;
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

        // Set convergence check to 0 (false) before the load steps start.
        stepConverged = 0;

        // Start load step
        for (int l = 0; l < (maxLoadSteps + 1); l++) // The +1 is so that the system is for conergence check and exiting later. 
        {
            printf("\n\tStep: %d, Load Step: %d . . . Running . . . \n\n", k, l);
            stepCounter += 1;

            // Init / reset global stiffness matrix and global internal force vetor for the current load step
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

            // Init / Reset increment displacement vector.
            for (int i = 0; i < Km; i++)
            {
                *(du + i) = 0;
            }

            // Loop through all elements
            for (int e = 0; e < nelements; e++)
            {
                printf("\t\tStep: %d, Load Step: %d . . . Analysing element #%d . . . \n", k, l, e);

                // Init / reset K_e and F_int_e for the current element
                initElementStiffnessMatrix(Ke, Kem);
                initElementInternalForceVector(FintE, Kem);

                // Get the displacements for the element's nodes, u_e, and store them. One value for each DOF.
                printf("\t\t\tDisplacements, u_e%d, at (%p):   ", e, (void*)elements[e].nodeids);
                for (int i = 0; i < nnodesElement; i++)
                {
                    if (elements[e].nodeids == NULL) {
                        printf("Error: elements[%d].nodeids is NULL\n", e);
                        exit(1);
                    }
                    int node_id = elements[e].nodeids[i];
                    if (node_id < 0 || node_id >= nnodes) {
                        printf("Error: node_id %d out of bounds for element %d (i=%d)\n", node_id, e, i);
                        exit(1);
                    }

                    for (int d = 0; d < DOF; d++) // Loop over each DOF
                    {
                        *(ue + i * DOF + d) = *(u + elements[e].nodeids[i] * DOF + d); // Store the displacement from the node with this id.
                    }
                }
                printPreciseMatrix(ue, 1, Kem);

                // Begin looping over the Gauss Points of the elements, each loop is for the i-th Gauss Point.
                for (int i = 0; i < elements->gp; i++)
                {
                    printf("\t\t\tElement #%d . . . Gauss Point #%d \n", e, i);

                    // Find the Jacobian, B matrix, and elastic D_e matrix for this Gauss Point.
                    quadJacobian(J, Jinv, Jn, detJ, sf[i].Ni, sf[i].NiPxi, sf[i].NiPeta, elements[e].coords, DOF);
                    quadBMatrix(B, Btrans, Jinv, sf[i].NiPxi, sf[i].NiPeta);
                    elasticDMatrixPlaneStrain(D, E, v);

                    // Get the strain for this Gauss Point with the displacement and B matrix.
                    displacementStrain(epsilon, B, ue, nnodesElement, DOF);

                    printf("\t\t\t\tStrain: \t");
                    for (int i = 0; i < Dn; i++)
                    {
                        printf("%.7f\t", *(epsilon + i));
                    }
                    printf("\n");

                    // Make trial stress
                    trialStress(sigmaTrial, D, epsilon, elements[e].epsilonP, Dn, i);

                    printf("\t\t\t\tTrial stress: ");
                    for (int i = 0; i < Dn; i++)
                    {
                        printf("%.2f\t", *(sigmaTrial + i));
                    }
                    printf("\n");

                    // Find deviatoric stress, von Mises equivalent stress, yield stress adjusted for existing plastic strain, and find von Mises Yield Function
                    deviatoricStress2D(sDeviatoric, sigmaTrial);
                    sigmaEq = vonMisesEquivalentStress2D(sDeviatoric);
                    sigmaYield = plasticCorrectedYieldStress(sigmaYieldInitial, H, elements[e].epsilonBarP[i]);

                    // Check to see if sigmaEq or sigmaYield are inf or are nan. If so, print error and exit the program.
                    if (isinf(sigmaEq) || isnan(sigmaEq) || isinf(sigmaYield) || isnan(sigmaEq))
                    {
                        fprintf(stderr, "Error: sigmaEq and/or sigmaYield are invalid. (inf/NaN). \n ");
                        exit(1);
                    }

                    // Check von Mises Yield Function for yielding of the node. If elastic, use elastic material stiffness matrix. If plastic, run return mapping.
                    f = vonMisesYieldFunction(sigmaEq, sigmaYield);

                    if (f > 0) // f > 0 -> yield
                    {
                        // Return mapping processes
                        printf("\t\t\t\t\tYIELD at GP:%d, \tEq Stress: %.4f, \tCorrected Yield Stress: %.4f\n", i, sigmaEq, sigmaYield);

                        // Unit deviatoric stress, n
                        unitDeviatoricStress(nUnitDeviatoric, sDeviatoric, Dn);
                        printf("\t\t\t\t\tn - unit deviatoric stress: ");
                        for (int i = 0; i < Dn; i++)
                        {
                            printf("%.2f\t", *(nUnitDeviatoric + i));
                        }
                        printf("\n");

                        // Elasto-plastic material stiffness matrix, D_ep
                        elastoPlasticDMatrix(D, nUnitDeviatoric, H, Dn);

                        // Plastic multiplier
                        deltaGamma = plasticMultiplier(f, G, H);

                        // Corrected trial stress
                        plasticStressCorrection(sigma, sigmaTrial, nUnitDeviatoric, G, deltaGamma, Dn);

                        printf("\t\t\t\t\tCoorected stress: ");
                        for (int i = 0; i < Dn; i++)
                        {
                            printf("%.2f\t", *(sigma + i));
                        }
                        printf("\n");

                        printf("\t\t\t\t\tG - shear modulus: %.4f   f - yield function %.4f   deltaGamma: %.10f \n", G, f, deltaGamma);

                        // Append the contribution to the element internal force vector, f^e_int, with the corrected stress. The element internal force vector is summed for each element over all the Gauss Points.
                        elementInternalForceVector(FintE, Btrans, sigma, *detJ, sf[i].weight, Kem, Dn);

                        // Trial plastic strain tensor
                        trialPlasticStrain(trialEpsilonP, elements[e].epsilonP, nUnitDeviatoric, deltaGamma, Dn, i);

                        // Trial plastic equivalent strain
                        elements[e].trialEpsilonBarP[i] = elements[e].epsilonBarP[i] + trialEquivalentPlasticStrain(deltaGamma);

                        // Store trial values in element, but they are not to be commited before the load step converges.
                        for (int j = 0; j < Dn; j++)
                        {
                            elements[e].trialSigma[i * elements[e].gp + j] = *(sigma + j);
                            elements[e].trialEpsilonP[i * elements[e].gp + j] = *(trialEpsilonP + j);
                        }
                        
                        
                        printf("\n\t\t\t\ttrialSigma: ");
                        for (int j = 0; j < Dn; j++) printf("%.1f\t", elements[e].trialSigma[i * elements[e].gp + j]);
                        printf("\n\t\t\t\tcommitedSigma: ");
                        for (int j = 0; j < Dn; j++) printf("%.1f\t", elements[e].sigma[i * elements[e].gp + j]);
                        printf("\n\t\t\t\ttrialEpsilonP: ");
                        for (int j = 0; j < Dn; j++) printf("%.6f\t", elements[e].trialEpsilonP[i * elements[e].gp + j]);
                        printf("\n\t\t\t\tcommitedEpsilonP: ");
                        for (int j = 0; j < Dn; j++) printf("%.6f\t", elements[e].epsilonP[i * elements[e].gp + j]);
                        printf("\n\t\t\t\ttrialEpsilonBarP: %.6f", elements[e].trialEpsilonBarP[i]);
                        printf("\n\t\t\t\tcommittedEpsilonBarP: %.6f\n", elements[e].epsilonBarP[i]);
                    }
                    else
                    {
                        // Append the contribution to the element internal force vector, f^e_int, with the trial stress. The element internal force vector is summed for each element over all the Gauss Points.
                        printf("\t\t\t\t\tNO Yield at GP:%d, \tEq Stress: %.4f, \tCorrected Yield Stress: %.4f\n", i, sigmaEq, sigmaYield);
                        elementInternalForceVector(FintE, Btrans, sigmaTrial, *detJ, sf[i].weight, Kem, Dn);
                    }

                    // Calculate and append the contribution to the element stiffness matrix, Ke, for the current Gauss Point. The element stiffness matrix is summed for each element over all the Gauss Points.
                    elementStiffnessMatrix(Ke, nnodesElement, DOF, B, Btrans, D, *detJ, sf[i].weight);
                }

                // Assemble into global stiffness matrix.
                if (skylineSolver)
                {
                    // Assemble into skyline matrix
                    globalSkylineStiffnessMatrix(Kskyline, Ke, elements[e].nodeids, DOF, nnodesElement, Kem);
                    applyFixedDisplacementSkylineStiffnessMatrix(uFixed, Kskyline);
                }
                {
                    // Assemble into regular matrix
                    globalStiffnessMatrix(K, Ke, elements[e].nodeids, DOF, nnodesElement, Km, Kem);
                    applyFixedDisplacementStiffnessMatrix(uFixed, K, Km);
                }
                
                // Assemble into global internal force vector.
                globalInternalForceVector(Fint, FintE, elements[e].nodeids, DOF, nnodesElement, Km, Kem);
            }

            // Print internal force vector
            printf("F_int - Internal force vector\n");
            for (int i = 0; i < Km; i++)
            {
                if (*(Fint + i) < 0) printf("%.4f   ", *(Fint + i));
                else printf("+%.4f   ", *(Fint + i));
            }
            printf("\n\n");

            // Compute the residual, r, as f_ext - f_int, and print f_ext_step
            printf("F_step - External force vector at the current step\n");
            for (int i = 0; i < Km; i++)
            {
                printf("%.4f   ", *(Fstep + i));
                *(r + i) = *(Fstep + i) - *(Fint + i);  
            }
            printf("\n\n");

            // Set fixed displacements in the residual vector
            applyFixedDisplacementResidualVector(uFixed, r, Km);

            // Print the residual vector
            printf("r - residual force vector at the current step\n");
            for (int i = 0; i < Km; i++)
            { 
                printf("%.4f   ", *(r + i));
            }
            printf("\n\n");

            // Check for convergence
            stepConverged = 1; // Make this prediction now, and correct it if ANY of the values in the residual are not below the threshold.
            for (int i = 0; i < Km; i++)
            {
                if ( fabs( *(r + i) ) > residualThreshold )
                {
                    stepConverged = 0; // If the current residual value is greater than the threshold, the solution is not converged.
                }
            }

            // If not converged, solve the system.
            if (!stepConverged)
            {
                // Check if this is the last load step, if so and the system did not converge, exit the program.
                if (l == maxLoadSteps)
                {
                    fprintf(stderr, "Error: System did not converge for Step# %d, Load Step# %d", k, l);
                    exit(1);
                }

                printf("\tStep: %d, Load Step: %d . . . NOT CONVERGED . . . SOLVING SYSTEM\n", k, l);

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

                printf("\t\tSolution of the system (transposed for ease of reading): \n\t\t");
                for (int i = 0; i < Km; i++) printf(" %.12f   ", *(du + i));
                printf("\n");

                // Add the displacement increment to the global displacement vector
                printf("\t\tTotal displacement of the system (transposed for ease of reading): \n\t\t");
                for (int i = 0; i < Km; i++)
                {
                    *(u + i) += *(du + i);
                    printf(" %.12f   ", *(u + i));
                }
                printf("\n\n");
            }
            else
            {
                commitTrialValuesAtGaussPoints(elements, nelements); // Commit trial values to commited values.
                printf("\tStep: %d, Load Step: %d . . . CONVERGED . . . \n\n", k, l);
                break;
            }
        }
    }

    printf("\n\nNewton-Raphson iteration has concluded after %d steps! \n\n", stepCounter);

    

    // Post-processing
    if (solutionFound)
    {
        //postProcessingElastoPlastic (u, element, sf, Km);
    }
    else
    {
        // No solution was found. Print error message.
        printf("No solutions were found for the analysis.");
    }


    // Freeing some variables from memory now that they are no longer needed.
    free(elements); 
    free(sf); 
    free(J); 
    free(detJ); 
    free(Jinv); 
    free(B); 
    free(Btrans);
    free(D);
    free(u);
    free(du);
    free(ue);
    free(Ke);
    free(K);
    free(Kskyline);
    free(r);
    free(Fint);
    free(Fext);
    free(FintE);
    free(epsilon);
    free(sigma);
    free(sigmaTrial);
    free(sDeviatoric);
    free(nUnitDeviatoric);
   
    return 0;
}