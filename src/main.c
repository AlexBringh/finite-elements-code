//
//	Author: Alexander B. Ringheim
//	Project: FEM in Plasticity
//	Date of creation: 16.10.2024
//

#include <stdio.h>

#include "flowRule.h"
#include "workHardening.h"
#include "skyline.h"
#include "solveCroutSkyline.h"

int main (char *args)
{
    /*
        The entry point of the program.
    */
    printf("Running Finite Elements Method! \n");

    // 1: Initialize -> Material Properties -> Stress and Strain Variables -> Initial Conditions and Read Input Data (Geometry, Boundary Conditions & Loading)

    // 2: Elastic Predictor -> Calculate Trial Stress (Assuming All Strains Are Initially Elastic)

    // 3: Yield Check  -> Evaluate the Yield Function Using The Trial Stress. -> Check if the Yield Condition is met.

    // 4: Plastic Corrector -> If the Yield Condition is met, apply the Flow Rule. Update the stress using Plastic Strain Increments. Compute the Corrector Step which involves Iterative Computation to satisfy the Yield Criterion. -> Integrate the Hardening Model. Update the Yield Surface for Isotropic Hardening or translate it for Kinematic Hardening.

    // 5: Stiffness Matrix Assembly. -> Form the Global Stiffness Matrix considering updated material properties. -> Include effects of plasticity in the Local Stiffness Matrix calculations.

    // 6: Solve System of Equations. -> Apply Boundary Conditions. -> Solve the Global System to find displacement increments.

    // 7: Post-Processing. -> Update Strains and Stresses from Displacement Increments. -> Check for Convergence and iterate if needed. -> Extract results such as Deformation, Stress Distributions, etc.

    return 0;
}