#include <stdio.h>
#include <stdlib.h>

#include "testCroutReduction.h"
#include "main.h"
#include "matrixUtils.h"

void testCroutReduction ();

int main ()
{
    printf("\nRUNNING TESTS FOR: Crout Reduction. \n");
    printf("croutReduction.c . . . \n");

    testCroutReduction();

    return 0;
}

void testCroutReduction ()
{
    int o = 2;
    double *A = malloc(o * o * sizeof(double));
    if (A == NULL) 
    {
        fprintf(stderr, "Error: Memory allocation for array A failed!");
        exit(1);
    }

    double *B = malloc(o * sizeof(double));
    if (B == NULL)
    {
        fprintf(stderr, "Error: Memory allocation for array B failed!");
        exit(1);
    }

    double *x = malloc(o * sizeof(double));
    if (x == NULL)
    {
        fprintf(stderr, "Error: Memory allocation for array x failed!");
        exit(1);
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
        exit(1);
    }

    B = malloc(o * sizeof(double));
    if (B == NULL)
    {
        fprintf(stderr, "Error: Memory allocation for array B failed!");
        exit(1);
    }

    x = malloc(o * sizeof(double));
    if (x == NULL)
    {
        fprintf(stderr, "Error: Memory allocation for array x failed!");
        exit(1);
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

    o = 3;
    A = malloc(o * o * sizeof(double));
    if (A == NULL) 
    {
        fprintf(stderr, "Error: Memory allocation for array A failed!");
        exit(1);
    }

    B = malloc(o * sizeof(double));
    if (B == NULL)
    {
        fprintf(stderr, "Error: Memory allocation for array B failed!");
        exit(1);
    }

    x = malloc(o * sizeof(double));
    if (x == NULL)
    {
        fprintf(stderr, "Error: Memory allocation for array x failed!");
        exit(1);
    }

    A[0] = 1;
    A[1] = 1;
    A[2] = 1;

    A[3] = 0;
    A[4] = 2;
    A[5] = 5;

    A[6] = 2;
    A[7] = 5;
    A[8] = -1;

    B[0] = 6;
    B[1] = -4;
    B[2] = 27;

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
        exit(1);
    }

    B = malloc(o * sizeof(double));
    if (B == NULL)
    {
        fprintf(stderr, "Error: Memory allocation for array B failed!");
        exit(1);
    }

    x = malloc(o * sizeof(double));
    if (x == NULL)
    {
        fprintf(stderr, "Error: Memory allocation for array x failed!");
        exit(1);
    }

    A[0] = 2;
    A[1] = 4;
    A[2] = -2;

    A[3] = 4;
    A[4] = 9;
    A[5] = -3;

    A[6] = -2;
    A[7] = -3;
    A[8] = 7;

    B[0] = 2;
    B[1] = 8;
    B[2] = 10;

    croutReduction(A, o, x, B);

    free(A);
    free(B);
    free(x);

    A = NULL;
    B = NULL;
    x = NULL;
}