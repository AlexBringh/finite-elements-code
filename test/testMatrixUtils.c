#include "testMatrixUtils.h"
#include <stdio.h>
#include <stdlib.h>

void testPrintMatrix();

int main ()
{
    // Test printMatrix()
    printf("\nRUNNING TESTS FOR: Matrix Utils: printMatrix(). \n");
    printf("matrixUtils.c . . . \n");
    testPrintMatrix();

    // Test augmentedMatrix()
    printf("\nRUNNING TESTS FOR: Matrix Utils: augmentedMatrix(). \n");
    printf("matrixUtils.c . . . \n");
    testAugmentedMatrix();

    return 0;
}

void testPrintMatrix ()
{
    // Make matrix A
    int m = 3;
    int n = m;
    double *A = malloc(m * n * sizeof(double));
    if (A == NULL)
    {
        fprintf(stderr, "Error (testAugmentedMatrix): Memory allocation for array 'A' failed! \n");
        exit(1);
    }

    A[0] = 1;
    A[1] = 2;
    A[2] = 3;
    A[3] = 4;
    A[4] = 5;
    A[5] = 6;
    A[6] = 7;
    A[7] = 8;
    A[8] = 9;

    // Test printMatrix()
    printMatrix(A, m, n);
}

void testAugmentedMatrix ()
{
    
    int m = 3;
    int n = m + 1;

    // Make augmented matrix Ab
    double *Ab = malloc(m * n * sizeof(double));
    if (Ab == NULL)
    {
        fprintf(stderr, "Error (testAugmentedMatrix): Memory allocation for array 'Ab' failed! \n");
        exit(1);
    }
    // Initialize memory slots of Ab to be 0.
    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < n; j++)
        {
            *(Ab + i*n + j) = 0;
        }
    }

    // Make matrix A
    double *A = malloc(m * m * sizeof(double));
    if (A == NULL)
    {
        fprintf(stderr, "Error (testAugmentedMatrix): Memory allocation for array 'A' failed! \n");
        exit(1);
    }

    A[0] = 1;
    A[1] = 2;
    A[2] = 3;
    A[3] = 4;
    A[4] = 5;
    A[5] = 6;
    A[6] = 7;
    A[7] = 8;
    A[8] = 9;

    // Make vector B
    double *B = malloc(m * sizeof(double));
    if (B == NULL)
    {
        fprintf(stderr, "Error (testAugmentedMatrix): Memory allocation for array 'B' failed! \n");
        exit(1);
    }

    B[0] = 200;
    B[1] = 500;
    B[2] = 800;

    // Test augmentedMatrix()
    augmentedMatrix(Ab, A, B, m);

    // Print input matrices
    printMatrix(A, m, m);
    printMatrix(B, m, 1);
    // Print result
    printMatrix(Ab, m, n);
}