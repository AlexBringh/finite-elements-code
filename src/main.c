//
//	Author: Alexander B. Ringheim
//	Project: FEM in Plasticity
//	Date of creation: 16.10.2024
//

#include <stdio.h>

#include "elements.h"
#include "jacobian.h"
#include "stiffnessMatrix.h"
#include "flowRule.h"
#include "workHardening.h"
#include "skyline.h"
#include "solveCroutSkyline.h"

int main (char *args)
{
    /*
        The entry point of the program.
    */
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
            - if f > 0 -> Compute delta gamma = f_yield / (3 G + H)
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

    // Initialize -> Material Properties -> Stress and Strain Variables -> Initial Conditions and Read Input Data (Geometry, Boundary Conditions & Loading)

    
    // Constant values for the analysis
    int gp = 4; // Gauss Points per element
    int nelements = 100; // Number of elements in the mesh.
    int nodeselement = 4; // Nodes per element
    int nnodes = 200; // Number of nodes in the mesh.
    int DOF = 2; // Degrees of freedom per node
    int Km = nnodes * DOF; // Number of rows / columns of the global stiffness matrix.

    // Allocate memory for the 'element' structs.
    quadElement *element = malloc(nelements * sizeof(quadElement));

    // Allocate memory for the global-, increment- and element displacement vector.
    double *u = malloc(nnodes * DOF * sizeof(double));
    double *du = malloc(nnodes * DOF * sizeof(double));
    double *ue = malloc(nodeselement * DOF * sizeof(double));

    // Get the shape functions for 2D Quad elements
    int j_m = 2; // Representing the size, m, in a m*m Jacobian matrix.
    double *Ni = malloc(j_m * j_m * sizeof(double));
    double *NiPxi = malloc(j_m * j_m * sizeof(double));
    double *NiPeta = malloc(j_m * j_m * sizeof(double));
    double *weights = malloc(gp * sizeof(double));
    quadShapeFunction(Ni, NiPxi, NiPeta, weights); // Gets the shape functions for a 2D quad element. Needs only be ran once.

    // Allocate memory for the Jacobian matrix, its determinant and the inverse.
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
    int Kem = gp * DOF;
    int Ken = Kem;
    double *Ke = malloc(Kem * Ken * sizeof(double));

    // Check how many nodes and degrees of freedom there are. If less than threshold, use normal matrix. If over threshold, use Skyline.
    // Allocate memory for the global stiffness matrix, (using Skyline).


    // Variables used in the Newton-Raphson iteration.
    int nodeID;

    // Newton-Raphson Iterator
    printf("Starting Newton-Raphson iteration.");
    int maxIterations = 100;
    for (int k = 0; k < maxIterations; k++)
    {
        printf("Running step: %1d . . . \n", k);

        // Apply load increment

        // Loop through all elements
        for (int e = 0; e < nelements; e++)
        {
            printf("\tIn step: %1d . . . Analysing element #%1d . . . \n", k, e);

            // Init K_e for current element
            for (int i = 0; i < Kem; i++)
            {
                for (int j = 0; j < Kem; j++)
                {
                    *(Ke + i * Kem + j) = 0; // Sets all the values of the element stiffness matrix to 0.
                }
            }

            // Get the displacements for the element's nodes, u_e, and store them. One value for each DOF.
            for (int i = 0; i < nodeselement * DOF; i++)
            {
                nodeID = *(element->nodeids + i); // Get the i-th node id stored in the element.
                *(ue + i) = *(u + nodeID); // Store the displacement from the node with this id.
            }

            // Begin looping over the Gauss Points of the elements
            for (int i = 0; i < element->gp; i++)
            {
                printf("\t\tElement #%1d . . . Gauss Point #%1d", e, i);


            }

        }
    }

    for (int i = 0; i < nelements; i++) // Loop through each element and find the local Jacobian matrix, its determinant, the B-matrix, the B-transpose, the D-matrix, and calculate the element stiffness matrix. 
    {
        double *xi; // TODO: Make logic for extracting nodes coordinates into xi and yi for each element.
        double *yi;
        quadJacobian(J, Jinv, j_m, detJ, Ni, NiPxi, NiPeta, xi, yi);
        quadBMatrix (B, Btrans, Jinv, NiPxi, NiPeta);
        // QUADELASTICDMATRIXPLANESTRAIN
        // QUADELASTOPLASTICDMATRIX
        quadElementStiffnessMatrix(Ke, gp, DOF, B, Btrans, D, *detJ, weights);

    }

    // Post-processing

    // Freeing some variables from memory now that they are no longer needed.
    free(element); // Element structs
    free(Ni); // Local shape functions.
    free(NiPxi); // Partial derivatives of local shape functions, of Xi (greek letter).
    free(NiPeta); // Åarial derivatives of local shape functions, of Eta.
    free(J); // Jacobian matrix.
    free(detJ); // Determinant of the Jacobian matrix.
    free(Jinv); // Inverse matrix of the Jacobian.
    free(B); // B matrix
    free(Btrans); // B transpose
    free(D);
    free(Ke); // Ke. element stiffness matrix

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
        - If f > 0 -> Plasticity (comput delta gamma, n, and D_ep), if f <= 0 -> elasticity (use D_e)
            - if f > 0 -> Compute delta gamma = f_yield / (3 G + H)
            - G = E / (2(1+v))
            - Correct stress on Gauss point as sigma = sigma_trial - 2 G * delta gamma * n
    
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
        - Run post-processing to extra results in useful formats.
    */

    return 0;
}