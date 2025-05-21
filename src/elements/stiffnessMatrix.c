#include <stdio.h>
#include <stdlib.h>

#include "stiffnessMatrix.h"
#include "matrixArithmetic.h"
#include "elements.h"


int quadBMatrix (double *B, double *Btrans, double *Jinv, double *NiPxi, double *NiPeta)
{
    /*
        Calculates the B matrix and B-transpose matrix from the values of J-inverse and the partial derivatives of the shape functions at the current Gauss Point.
        This function is ran once for each Gauss Point for each element.

        Input:
        double *B      -> Pointer to B matrix. (must be allocated memory for a 3x(2 * nnodesElement) matrix, layed out in array form). Results are stored here.
        double *Btrans -> Pointer to B-transpose matrix. (must be allocated memory for a (2 * nnodesElement)x3 matrix, layed out in array form). Results are stored here.
        double *Jinv   -> Pointer to the calculated J-inverse matrix for the current Gauss Point.
        double *NiPxi  -> Pointer to the partial derivatives with respects to xi of the shape functions.
        double *NiPeta -> Pointer to the partial derivatives with respects to eta of the shape functions.
    */
    int nnodes = 4;
    int dim = 2;

    double NiPx;
    double NiPy;

    // Init / reset B-matrix and B-transpose
    for (int i = 0; i < nnodes * dim * 3; i++)
    {
        *(B + i) = 0;
        *(Btrans + i) = 0;
    }


    // Calculate NiPx and NiPy for all 4 nodes of the element. Assemble the B-matrix and B-transpose for the Gauss Point.
    for (int i = 0; i < nnodes; i++)
    {
        NiPx = *(Jinv + 0) * *(NiPxi + i) + *(Jinv + 1) * *(NiPeta + i); // Calculate the partial derivative for the global x coordinate of the shape functions.
        NiPy = *(Jinv + 2) * *(NiPxi + i) + *(Jinv + 3) * *(NiPeta + i); // Calculate the partial derivative for the global y coordinate of the shape functions.

        *(B + i * 2) = NiPx; // First row of the B-matrix.
        *(B + 9 + i * 2) = NiPy; // Second row of the B-matrix.
        *(B + 16 + i * 2) = NiPy; // Third row, first value of the B-matrix.
        *(B + 17 + i * 2) = NiPx; // Third row, second value of the B-matrix.

        *(Btrans + i * 3 * 2) = NiPx; // First column of the B-transpose.
        *(Btrans + i * 3 * 2 + 4) = NiPy; // Second column of the B-transpose.
        *(Btrans + i * 3 * 2 + 2) = NiPy; // Third column, first value of the B-transpose.
        *(Btrans + i * 3 * 2 + 5) = NiPx; // Third column, second value of the B-transpose.
    }

    return 0;
}


int elasticDMatrixPlaneStress (double *D, double E, double v)
{
    /*
        Good for thin structures where the thickness is small compared to other dimensions and loads are in-plane.
        This function is used for systems where Plane Stress is assumed.

        The elastic material stiffness matrix is used to relate stress and strain by Hooke's law. 
        The elastic matrix is used for the domain where the material can be assumed not the undergo permanent deformation.

        Input:
        double *D -> Pointer to D matrix (must be allocated memory for a 3x3 matrix, layed out in array form). Results / assembled D_e matrix is stored here.
        double E  -> Young's Modulus
        double v  -> Poisson's Ratio
    */
   
    double multiplier = E  / (1 - v * v);

    *(D + 0) = multiplier * 1;
    *(D + 1) = multiplier * v;
    *(D + 2) = 0;
    *(D + 3) = multiplier * v;
    *(D + 4) = multiplier * 1;
    *(D + 5) = 0;
    *(D + 6) = 0;
    *(D + 7) = 0;
    *(D + 8) = multiplier * (1 - v) / 2;
}

int elasticDMatrixPlaneStrain (double *D, double E, double v)
{
    /*
        Good for thick- or long bodies where the deformation in one direction is negligible.
        This function is used for systems where Plane Strain is assumed.

        The elastic material stiffness matrix is used to relate stress and strain by Hooke's law. 
        The elastic matrix is used for the domain where the material can be assumed not the undergo permanent deformation.

        Input:
        double *D -> Pointer to D matrix (must be allocated memory for a 3x3 matrix, layed out in array form). Results / assembled D_e matrix is stored here.
        double E  -> Young's Modulus
        double v  -> Poisson's Ratio
    */
    double multiplier = E * (1 - v) / ((1 + v) * (1 - 2 * v));

    *(D + 0) = multiplier * 1;
    *(D + 1) = multiplier * v / (1 - v);
    *(D + 2) = 0;
    *(D + 3) = multiplier * v / (1 - v);
    *(D + 4) = multiplier * 1;
    *(D + 5) = 0;
    *(D + 6) = 0;
    *(D + 7) = 0;
    *(D + 8) = multiplier * (1 - 2 * v) / (2 * (1 - v));
}


int elementStiffnessMatrix (double *Ke, int nnodesElement, int DOF, double *B, double *Btrans, double *D, double detJ, double weight)
{
    /*
        Calculates the element stiffness matrix contribution for a given Gauss Point.
        This function will be ran once for every Gauss Point of an element, with its own unique B matrix, Jacobian, D matrix, and weight (though integration weight is not always unique. See 2D quad elements.)
        The function adds teh contribution to the element stiffness matrix for the current element. The matrix should be reset at the start of every element in the system, but not for every Gauss Point.

        Input: 
        double *Ke -> Element stiffness matrix
    */

    int m = nnodesElement * DOF;

    double *Ktemp1 = malloc (m * 3 * sizeof(double)); // Temp matrix(mx3) for Btrans(mx3) x D(3x3)
    double *Ktemp2 = malloc (m * m * sizeof(double)); // Temp matrix(mxm) for Ktemp1(mx3) x B(3xm)

    for (int j = 0; j < m * 3; j++) *(Ktemp1 + j) = 0; // Initialize / empty Ktemp1(mx3)
    for (int j = 0; j < m * m; j++) *(Ktemp2 + j) = 0; // Initialize / empty Ktemp2(mxm)

    matrixMultiply(Btrans, D, Ktemp1, 8, 3, 3);
    matrixMultiply(Ktemp1, B, Ktemp2, 8, 3, 8);
    
    // Summation of the element stiffness matrix for the current Gauss Point.
    for (int i = 0; i < m * m; i++)
    {
        // Multiply the det(J) and weight to each cell of Ktemp2. Add each cell of the current matrix to the element stiffness matrix, Ke.
        *(Ke + i) += *(Ktemp2 + i) * detJ * weight;
    }

    free(Ktemp1);
    free(Ktemp2);
}


int initElementStiffnessMatrix (double *Ke, int Kem)
{
    /*
        Initializes / resets all values of the Ke matrix to 0. Determines the size by the number of nodes per element, times the degrees of freedom per node.
    */
   
    for (int i = 0; i < Kem * Kem; i++)
    {
        *(Ke + i) = 0;
    }
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
    for (int i = 0; i < Km * Km; i++)
    {
        *(K + i) = 0;
    }
}


int globalStiffnessMatrix (double *K, double *Ke, int *nodeids, int dof, int nnodesElement, int Km, int Kem)
{
    /*
        Assembles the element stiffnes matrix into the global stiffness matrix. 
        Determines the correct indices based on the node id found in the current element.

        Inputs:
        double *K          ->  Pointer to global stiffness matrix. Results are stored here.
        double *Ke         ->  Pointer to element stiffness matrix. 
        int *nodeids       ->  Pointer to nodeids for the nodes building the current element.
        int dof            ->  Degrees of freedom of the nodes.
        int nnodesElement  ->  Number of nodes in the element.
        int Km             ->  Size, m_global, of the rows/columns of the global stiffness matrix.
        int Kem            ->  Size, m_element, of the rows/columns of the element stiffness matrix.
    */

    // Allocate memory for array to hold global degrees of freedom indices.
    int *globDOFs = malloc (nnodesElement * dof * sizeof(int));

    // Calculate and store global degrees of freedom indices from each of the nodes ids stored in the 'nodeids' array.
    for (int i = 0; i < nnodesElement; i++)
    {
        for (int d = 0; d < dof; d++)
        {
            *(globDOFs + i * dof + d) = *(nodeids + i) * dof + d;
        }
    }

    int Ki; // Row indices for K
    int Kj; // Col indices for K

    // Assemble the global stiffness matrix by looping over each element cell in the element stiffness matrix, 
    // finding the global indices corresponding, and adding the element cell values to the global cell values.
    for (int i = 0; i < Kem; i++)
    {
        Ki = *(globDOFs + i); // Row indices for K
        for (int j = 0; j < Kem; j++)
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

