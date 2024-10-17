//
//	Author: Alexander B. Ringheim
//	Project: FEM in Plasticity
//	Date of creation: 16.10.2024
//

#include <stdio.h>
#include <stdlib.h>
#include "main.h"
#include "gauss_elimination.h"
#include "matrix_utils.h"

int main (char *args)
{
    printf("Hello World! \n");

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

    double *C = (double *) malloc(m * h * sizeof(double));
    if (C == NULL)
    {
        fprintf(stderr, "Memory allocation for array C failed! \n");
        return -1;
    }


    int i = 3;
    int j = 4;
    double *D = (double *) malloc(i * j * sizeof(double));
    if (D == NULL)
    {
        fprintf(stderr, "Memory allocation for array D failed! \n");
        return -1;
    }

    D[0] = 1;
    D[1] = 1;
    D[2] = 1;
    D[3] = 6;

    D[4] = 0;
    D[5] = 2;
    D[6] = 5;
    D[7] = -4;

    D[8] = 2;
    D[9] = 5;
    D[10] = -1;
    D[11] = 27;

    double *E = (double *) malloc(i * sizeof(double));
    if (E == NULL)
    {
        fprintf(stderr, "Memory allocation for array E failed! \n");
        return -1;
    }
    for (int it = 0; it < j; it++)
    {
        E[it] = 0;
    }

    gauss_elimination(D, i, E);


    // Clear variables from memory.
    free(A);
    free(B);
    free(C);
    free(D);
    free(E);

    // Set variables to NULL.
    A = NULL;
    B = NULL;
    C = NULL;
    D = NULL;
    E = NULL;

    return 0;
}