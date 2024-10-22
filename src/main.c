//
//	Author: Alexander B. Ringheim
//	Project: FEM in Plasticity
//	Date of creation: 16.10.2024
//

#include <stdio.h>
#include <stdlib.h>
#include "main.h"
#include "gaussElimination.h"
#include "matrixUtils.h"

int main (char *args)
{
    printf("\n");

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