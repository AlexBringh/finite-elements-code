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
    printf("\n");

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

    


    int m = 3; // Rows for D-matrix and E-vector.
    int n = 4; // Columns for D-matrix
    double *D = malloc(m * n * sizeof(double));
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

    double *E = malloc(m * sizeof(double));
    if (E == NULL)
    {
        fprintf(stderr, "Memory allocation for array E failed! \n");
        return -1;
    }
    for (int i = 0; i < m; i++)
    {
        E[i] = 0;
    }

    gaussElimination(D, m, E);

    // Free variables after use
    free(D);

    m = 3; // Rows for F-matrix and G-vector.
    n = 4; // Columns for F-matrix
    double *F = malloc(m * n * sizeof(double));
    if (F == NULL)
    {
        fprintf(stderr, "Memory allocation for array F failed! \n");
        return -1;
    }
    
    // Make init variable and 
    F[0] = 2;
    F[1] = 4;
    F[2] = -2;
    F[3] = 2;

    F[4] = 4;
    F[5] = 9;
    F[6] = -3;
    F[7] = 8;
    
    F[8] = -2;
    F[9] = -3;
    F[10] = 7;
    F[11] = 10;

    double *G = malloc(m * sizeof(double));
    if (G == NULL)
    {
        fprintf(stderr, "Memory allocation for array E failed! \n");
        return -1;
    }
    for (int i = 0; i < m; i++)
    {
        G[i] = 0;
    }

    gaussElimination(F, m, G);

    free(F);

    // Free solution variables after use.
    free(E);
    free(G);
    D = NULL;
    E = NULL;
    F = NULL;
    G = NULL;

    return 0;
}