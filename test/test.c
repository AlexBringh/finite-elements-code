//
//	Author: Alexander B. Ringheim
//	Project: Test
//	Date of creation: 14.10.2024
//

#include "test.h"
#include "matrix_arithmetic.h"
#include "print_matrix.h"

int main (char* args)
{

    int m = 4;
    int n = 3;
    int h = 2;

    double *A = (double *) malloc(m * n * sizeof(double));
    if (A == NULL)
    {
        fprintf (stderr, "Memory allocation for array A failed! \n");
        return -1;
    }
    // Set values for the A-array.
    for (int i = 0; i < m * n; i++)
    {
        A[i] = i + 1;
    }


    double *B = (double *) malloc(n * h * sizeof(double));
    if (B == NULL)
    {
        fprintf (stderr, "Memory allocation for array B failed! \n");
        return -1;
    }
    // Set values for the B-array.
    for (int i = 0; i < n * h; i++)
    {
        B[i] = i + 1;
    }
    
    /*
    double A = {
         1,  2,  3,
         4,  5,  6,
         7,  8,  9,
        10, 11, 12
    };
    double B[6] = {
        1, 2,
        3, 4,
        5, 6
    };
    */
    

    double *C = (double *) malloc(m * h * sizeof(double));
    if (C == NULL)
    {
        fpritnf(stderr, "Memory allocation for array C failed! \n");
    }

    matrix_multiply(A, B, C, m, n, h);
    print_matrix(A, m, n);
    print_matrix(B, n, h);
    print_matrix(C, m, h);

    // Clear variables from memory.
    free(A);
    free(B);
    free(C);

    // Set variables to NULL.
    A = NULL;
    B = NULL;
    C = NULL;

    return 0;
}
