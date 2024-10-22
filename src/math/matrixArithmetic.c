//
//	Author: Alexander B. Ringheim
//	Project: Finite Element Method for Plasticity analysis of materials and structures.
//	Date of creation: 14.10.2024
//

#include "matrixArithmetic.h"

int matrixSum (double *A, double *B, double *C, int m, int n)
{
    /*
        Sums the values of the A-matrix and the B-matrix, and stores the result in the given C-matrix.
        A and B must be equal m x n (row x column) matrices, or the program will cause errors.
    */

    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < n; j++)
        {
            *(C + i * n + j) = *(A + i*n + j) + *(A + i*n + j);
        }
    }

    return 0;
}

int matrixSubtract (double *A, double *B, double *C, int m, int n)
{
    /*
        Subtracts the values of the A-matrix using the B-matrix and stores the result in the given C-matrix.
        A and B must be equal m x n (row x column) matrices, or the program will cause errors.
    */

    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < n; j++)
        {
            *(C + i * n + j) = *(A + i*n + j) - *(A + i*n + j);
        }
    }

    return 0;
}

int matrixMultiply (double *A, double *B, double *C, int m, int n, int h)
{
    /*
        Multiplies 2 matrices together (A x B) and stores the result in the matrix C as a m x h (row x column) matrix.
        A must be a m x n matrix (row x column), whilst B must be a n x o (row x column) matrix.
    */

    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < h; j++)
        {
            for (int k = 0; k < n; k++)
            {
                *(C + i*h + j) += *(A + i*n + k) * *(B + k*h + j);
            }
        }
    }

    return 0;
}
