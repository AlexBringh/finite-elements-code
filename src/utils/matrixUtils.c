#include "matrixUtils.h"

void printMatrix (double *A, int m, int n)
{
    /*
        Prints the matrix of m x n (row x column) dimension.
    */
    printf("Matrix (%dx%d) \n", m, n);

    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < n; j++)
        {
            double a = *(A + n*i + j);
            printf("%.4f\t", a);
        }
        printf("\n");
    }
    printf("\n");
}


int augmentedMatrix (double *Ab, double *A, double *B, int m)
{
    /*
        Makes augmented matrix out of intputed m x m-matrix A, and m x 1-vector B. 
        The result is stored in the inputed Ab matrix, which must be a m x (m + 1)-matrix.
    */
    
    int n = m + 1;

    // Build the augmented matrix
    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < n; j++)
        {
            if (j < m) // If it is less than m (rows), then we want a value from the original A.
            {
                *(Ab + i * n + j) = *(A + i*m + j);
            }
            else // If it is equal or higher than m (will only ever be equal), then we want a value from B.
            {
                *(Ab + i * n + m) = *(B + i);
            }
        }
    }

    return 0;
}