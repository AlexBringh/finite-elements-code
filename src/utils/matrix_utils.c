#include "matrix_utils.h"

int print_matrix (double *A, int m, int n)
{
    /*
        Prints the matrix of m x n (row x column) dimension.
    */
    printf("Matrix (%dx%d) at location: %p \n", m, n, A);

    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < n; j++)
        {
            double a = *(A + n*i + j);
            printf("%.2f\t", a);
        }
        printf("\n");
    }
    printf("\n");
}