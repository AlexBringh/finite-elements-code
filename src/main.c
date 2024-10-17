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

    double *A = malloc(m * n * sizeof(double));
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


    double *B = malloc(n * h * sizeof(double));
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

    double *C = malloc(m * h * sizeof(double));
    if (C == NULL)
    {
        fprintf(stderr, "Memory allocation for array C failed! \n");
        return -1;
    }

    // Free variables after use.
    free(A);
    free(B);
    free(C);
    A = NULL;
    B = NULL;
    C = NULL;


    int i = 3; // Rows for D-matrix and E-vector.
    int j = 4; // Columns for D-matrix
    double *D = malloc(i * j * sizeof(double));
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

    double *E = malloc(i * sizeof(double));
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

    // Free variables after use
    free(D);
    free(E);
    D = NULL;
    E = NULL;


    i = 3; // Rows for F-matrix and G-vector.
    j = 4; // Columns for F-matrix
    double *F = malloc(i * j * sizeof(double));
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

    double *G = malloc(i * sizeof(double));
    if (G == NULL)
    {
        fprintf(stderr, "Memory allocation for array E failed! \n");
        return -1;
    }
    for (int it = 0; it < j; it++)
    {
        G[it] = 0;
    }

    gauss_elimination(F, i, G);

    double *H = malloc(i * j * sizeof(double));


    // Free variables after use.
    free(F);
    free(G);
    F = NULL;
    G = NULL;

    return 0;
}