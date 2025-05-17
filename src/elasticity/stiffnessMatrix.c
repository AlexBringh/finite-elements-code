#include <stdio.h>
#include <stdlib.h>

#include "stiffnessMatrix.h"
#include "matrixArithmetic.h"
#include "elements.h"

int quadBMatrix (double *B, double *Btrans, double *Jinv, double *NiPxi, double *NiPeta)
{
    /*
        TODO: Documentation
    */
    int nnodes = 4;
    double *NiPx = malloc(nnodes * sizeof(double)); // Allocate memory for partial deriv. of shape func., N_i, for global node coord., x, for 4 nodes (for quad mesh.)
    double *NiPy = malloc(nnodes * sizeof(double)); // Allocate memory for partial deriv. of shape func., N_i, for global node coord., y, for 4 nodes (for quad mesh.)

    // Calculate NiPx and NiPy for all 4 nodes of the element.
    for (int i = 0; i < nnodes; i++)
    {
        *(NiPx + i) = *(Jinv + 0) * *(NiPxi) + *(Jinv + 1) * *(NiPeta + i);
        *(NiPy + i) = *(Jinv + 2) * *(NiPxi) + *(Jinv + 3) * *(NiPeta + i);
    }

    // Assemble the B-matrix for the element.
    for (int i = 0; i < nnodes; i++)
    {
        *(B + i * 2) = *(NiPx + i); // First row of the B-matrix.
        *(B + 9 + i * 2) = *(NiPy + i); // Second row of the B-matrix.
        *(B + 16 + i * 2) = *(NiPy + i); // Third row, first value of the B-matrix.
        *(B + 17 + i * 2) = *(NiPx + i); // Third row, second value of the B-matrix.
    }

    // Assemble the B-transpose for the element.
    for (int i = 0; i < nnodes; i++)
    {
        *(Btrans + i * 3 * 2) = *(NiPx + i); // First column of the B-transpose.
        *(Btrans + i * 3 * 2 + 4) = *(NiPy + i); // Second column of the B-transpose.
        *(Btrans + i * 3 * 2 + 2) = *(NiPy + i); // Third column, first value of the B-transpose.
        *(Btrans + i * 3 * 2 + 5) = *(NiPx + i); // Third column, second value of the B-transpose.
    }

    free(NiPx);
    free(NiPy);

    return 0;
}

int elasticDMatrixPlaneStress ()
{
    /*
        TODO: Documentaiton
        Good for thin structures where the thickness is small compared to other dimensions and loads are in-plane.
    */
}

int elasticDMatrixPlaneStrain ()
{
    /*
        TODO: Documentation
        Good for thick- or long bodies where the deformation in one direction is negligible.
    */

}

int elastoPlasticDMatrix ()
{
    /*
        TODO: Documentation
    */
}

int quadElementStiffnessMatrix (double *Ke, int gp, int DOF, double *B, double *Btrans, double *D, double detJ, double *weights)
{
    /*
        TODO: Documentation
    */
    // Initialize / empty the element stiffness matrix
    for (int i = 0; i < gp * DOF; i++)
    {
        *(Ke + i) = 0;
    }

    double *Ktemp1 = malloc (8 * 3 * sizeof(double)); // Temp matrix(8x3) for Btrans(8x3) x D(3x3)
    double *Ktemp2 = malloc (8 * 8 * sizeof(double)); // Temp matrix(8x8) for Ktemp1(8x3) x B(3x8)
    
    // Summation of the element stiffness matrix for each Gauss Point, i
    for (int i = 0; i < gp; i++)
    {
        for (int j = 0; j < 8*3; j++) *(Ktemp1 + j) = 0; // Initialize / empty Ktemp1(8x3)
        for (int j = 0; j < 8*8; j++) *(Ktemp2 + j) = 0; // Initialize / empty Ktemp2(8x8)

        matrixMultiply(Btrans, D, Ktemp1, 8, 3, 3);
        matrixMultiply(Ktemp1, B, Ktemp2, 8, 3, 8);

        // Multiply the det(J) and weight to each cell of Ktemp2. Add each cell of the current matrix to the element stiffness matrix, Ke.
        for (int j = 0; j < gp * DOF; j++)
        {
            *(Ke + j) += *(Ktemp2 + j) * detJ * *(weights + i);
        }
    }

    free(Ktemp1);
    free(Ktemp2);
}

int getElementDOFs (int dof, int nelements)
{

}

int globalStiffnessMatrixSkyline ()
{

}

int initGlobalStiffnessMatrix (double *K, int Km)
{
    /*
        Sets all values of the global stiffness matrix to 0. Used for initializing and clearing the global stiffness matrix.
        Cannot be used for Skylie Matrix structures.

        Inputs:
        double *K -> pointer to global stiffness matrix in array-form.K
        int Km -> number of columns (rows) of the matrix.
    */
    for (int i = 0; i < Km; i++)
    {
        for (int j = 0; j < Km; j++)
        {
            *(K + i * Km + j) = 0;
        }
    }
}

int globalStiffnessMatrix (double *K, double *Ke, quadElement *element, int Km, int Kem)
{
    /*
        TODO: Documentation
    */
    //

    // Dereference values from the current element.
    int dof = element->dof;
    int nnodes = element->nnodes;

    // Allocate memory for array to hold global degrees of freedom indices.
    int *globDOFs = malloc (nnodes * dof * sizeof(int));

    // Calculate and store global degrees of freedom indices from each of the nodes ids stored in the connectivity array.
    for (int i = 0; i < nnodes; i++)
    {
        for (int j = 0; j < dof; j++)
        {
            *(globDOFs + i * dof + j) = *(element->connectivity + i) * dof + j;
        }
    }

    int Ki; // Row indices for K
    int Kj; // Col indices for K

    // Assemble the global stiffness matrix
    for (int i = 0; i < Km; i++)
    {
        Ki = *(globDOFs + i); // Row indices for K
        for (int j = 0; j < Km; j++)
        {

            // Check that indices are not out of bounds for the global stiffness matrix.
            if (Ki >= Km || Kj >= Km) 
            {
                fprintf(stderr, "Index out of bounds: Ki=%d, Kj=%d\n", Ki, Kj);
                return 1;
            }

            Kj = *(globDOFs + j); // Column indices for K

            // K[Ki][Kj] += Ke[i][j];
            *(K + Ki * Km + Kj) += *(Ke + i * Kem + j);
        }
    }

    free(globDOFs);
}

