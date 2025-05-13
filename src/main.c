//
//	Author: Alexander B. Ringheim
//	Project: FEM in Plasticity
//	Date of creation: 16.10.2024
//

#include <stdio.h>

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

    // 1: Initialize -> Material Properties -> Stress and Strain Variables -> Initial Conditions and Read Input Data (Geometry, Boundary Conditions & Loading)

    // 2: Elastic Predictor -> Calculate Trial Stress (Assuming All Strains Are Initially Elastic).

    // 3: Yield Check  -> Evaluate the Yield Function Using The Trial Stress. -> Check if the Yield Condition is met.

    // 4: Plastic Corrector -> If the Yield Condition is met, apply the Flow Rule. Update the stress using Plastic Strain Increments. Compute the Corrector Step which involves Iterative Computation to satisfy the Yield Criterion. -> Integrate the Hardening Model. Update the Yield Surface for Isotropic Hardening or translate it for Kinematic Hardening.

    // 5: Stiffness Matrix Assembly. -> Form the Global Stiffness Matrix considering updated material properties. -> Include effects of plasticity in the Local Stiffness Matrix calculations.
    
    // Get the shape functions for 2D Quad elements
    int j_m = 2; // Representing the size, m, in a m*m Jacobian matrix.
    double *Ni = malloc(j_m * j_m * sizeof(double));
    double *NiPxi = malloc(j_m * j_m * sizeof(double));
    double *NiPeta = malloc(j_m * j_m * sizeof(double));
    quadShapeFunction(Ni, NiPxi, NiPeta); // Gets the local shape functions for a quad element. Needs only be ran once.

    // Allocate memory for the Jacobian matrix, its determinant and the inverse.
    double *J = malloc(j_m * j_m * sizeof(double)); // This is initialized with values in the Jacobian function.
    double *detJ = malloc(sizeof(double)); // Store the determinant value of the Jacobian here.
    double *Jinv = malloc (j_m * j_m * sizeof(double)); // This is also initialized in the Jacobian function.

    // Allocate memory for the B matrix, B transpose, D and the element stiffness matrix.
    int nlocnode = 4; // Number of local nodes per element.
    int Bm = 3; // Rows of the B matrix.
    int Bn = nlocnode * 2; // Each node adds 2 columns to the B-matrix.
    double *B = malloc(Bn * Bm * sizeof(double)); 
    double *Btrans = malloc(Bn * Bm * sizeof(double));

    // Allocate memory for the global stiffness matrix, using Skyline.

    int nelements = 100; // PLACEHOLDER VALUE.
    for (int i = 0; i < nelements; i++) // Loop through each element and find the local Jacobian matrix, its determinant, the B-matrix, the B-transpose, the D-matrix, and calculate the element stiffness matrix. 
    {
        double *xi; // TODO: Make logic for extracting nodes coordinates into xi and yi for each element.
        double *yi;
        quadJacobian(J, Jinv, j_m, detJ, Ni, NiPxi, NiPeta, xi, yi);
        quadBMatrix (B, Btrans, Jinv, NiPxi, NiPeta);
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


    // ASSEMBLE THE GLOBAL STIFFNESS MATRIX AFTER EVERY ELEMENT STIFFNESS MATRIX HAS BEEN FOUND.

    // 6: Solve System of Equations. -> Apply Boundary Conditions. -> Solve the Global System to find displacement increments.

    // 7: Post-Processing. -> Update Strains and Stresses from Displacement Increments. -> Check for Convergence and iterate if needed. -> Extract results such as Deformation, Stress Distributions, etc.

    return 0;
}