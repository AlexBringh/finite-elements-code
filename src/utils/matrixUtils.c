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
            printf("%.3f\t", a);
        }
        printf("\n");
    }
    printf("\n");
}