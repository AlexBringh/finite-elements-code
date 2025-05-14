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
    printf("Running Finite Elements Simulation on Plasticity Equations! \n");

    // Initialize -> Material Properties -> Stress and Strain Variables -> Initial Conditions and Read Input Data (Geometry, Boundary Conditions & Loading)

    
    // Constant values for the analysis
    int gp = 4; // Gauss Points per element
    int DOF = 2; // Degrees of freedom per node / Gauss Point

    // Get the shape functions for 2D Quad elements
    int j_m = 2; // Representing the size, m, in a m*m Jacobian matrix.
    double *Ni = malloc(j_m * j_m * sizeof(double));
    double *NiPxi = malloc(j_m * j_m * sizeof(double));
    double *NiPeta = malloc(j_m * j_m * sizeof(double));
    double *weights = malloc(gp * sizeof(double));
    quadShapeFunction(Ni, NiPxi, NiPeta, weights); // Gets the local shape functions for a quad element. Needs only be ran once.

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
    int Ken = gp * DOF;
    double *Ke = malloc(Kem * Ken * sizeof(double));

    // Allocate memory for the global stiffness matrix, using Skyline.


    // Newton-Raphson Iterator
    printf("Starting Newton-Raphson iteration.");
    int maxIterations = 100;
    for (int k = 0; k < maxIterations; k++)
    {
        printf("Running step: %1d . . . \n", k);

    }

    int nelements = 100; // PLACEHOLDER VALUE.
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

    // Freeing some variables from memory now that they are no longer needed.
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
    */

    return 0;
}