//
//	Author: Alexander B. Ringheim
//	Project: FEM in Plasticity
//	Date of creation: 16.10.2024
//

#include <stdio.h>
#include <stdlib.h>
#include "main.h"
#include "gaussElimination.h"
#include "croutReduction.h"
#include "matrixUtils.h"

int main (char *args)
{
    int o = 2;
    double *A = malloc(o * o * sizeof(double));
    if (A == NULL) 
    {
        fprintf(stderr, "Error: Memory allocation for array A failed!");
        return -1;
    }

    double *B = malloc(o * sizeof(double));
    if (B == NULL)
    {
        fprintf(stderr, "Error: Memory allocation for array B failed!");
        return -1;
    }

    double *x = malloc(o * sizeof(double));
    if (x == NULL)
    {
        fprintf(stderr, "Error: Memory allocation for array x failed!");
        return -1;
    }

    A[0] = 6;
    A[1] = 3;
    A[2] = 4;
    A[3] = 3;

    B[0] = 4;
    B[1] = 2;

    croutReduction(A, o, x, B);


    free(A);
    free(B);
    free(x);
    A = NULL;
    B = NULL;
    x = NULL;

    o = 3;
    A = malloc(o * o * sizeof(double));
    if (A == NULL) 
    {
        fprintf(stderr, "Error: Memory allocation for array A failed!");
        return -1;
    }

    B = malloc(o * sizeof(double));
    if (B == NULL)
    {
        fprintf(stderr, "Error: Memory allocation for array B failed!");
        return -1;
    }

    x = malloc(o * sizeof(double));
    if (x == NULL)
    {
        fprintf(stderr, "Error: Memory allocation for array x failed!");
        return -1;
    }

    A[0] = 3;
    A[1] = -0.1;
    A[2] = -0.2;

    A[3] = 0.1;
    A[4] = 7;
    A[5] = -0.3;

    A[6] = 0.3;
    A[7] = -0.2;
    A[8] = 10;

    B[0] = 5;
    B[1] = -2;
    B[2] = 9;

    croutReduction(A, o, x, B);

    free(A);
    free(B);
    free(x);

    A = NULL;
    B = NULL;
    x = NULL;

    return 0;
}