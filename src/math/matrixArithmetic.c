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

        Inputs:
        double *A -> Pointer to A matrix (m x n)
        double *B -> Pointer to B matrix (m x n)
        double *C -> Pointer to C matrix (m x n) (results are stored here)
        int m -> size m
        int n -> size n
    */

    for (int i = 0; i < m * n; i++)
    {
        *(C + i) = *(A + i) + *(B + i);
    }

    return 0;
}

int matrixSubtract (double *A, double *B, double *C, int m, int n)
{
    /*
        Subtracts the values of the A-matrix using the B-matrix and stores the result in the given C-matrix.
        A and B must be equal m x n (row x column) matrices, or the program will cause errors.

        Inputs:
        double *A -> Pointer to A matrix (m x n)
        double *B -> Pointer to B matrix (m x n)
        double *C -> Pointer to C matrix (m x n) (results are stored here)
        int m -> size m
        int n -> size n
    */

    for (int i = 0; i < m * n; i++)
    {
        *(C + i) = *(A + i) - *(B + i);
    }

    return 0;
}

int matrixMultiply (double *A, double *B, double *C, int m, int n, int h)
{
    /*
        Multiplies 2 matrices together (A x B) and stores the result in the matrix C as a m x h (row x column) matrix.
        A must be a m x n matrix (row x column), whilst B must be a n x h (row x column) matrix.

        Inputs:
        double *A -> Pointer to A matrix (m x n)
        double *B -> Pointer to B matrix (n x h)
        double *C -> Pointer to C matrix (m x h) (results are stored here)
        int m -> size m
        int n -> size n
        int h -> size h
    */

    // Zero out C
    for (int i = 0; i < m * h; i++)
    {
        *(C + i) = 0;
    }

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
