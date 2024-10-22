//
//	Author: Alexander B. Ringheim
//	Project: Finite Element Method solver for Plasticity analysis of materials and structures.
//	Date of creation: 21.10.2024
//

#include "croutReduction.h"
#include "matrixUtils.h"

int croutReduction (double *A, int m, double *x)
{
    int n = m; // n represents the columns of the A-matrix.

    // Allocate memory slots for the L-matrix.
    double *L = malloc(m * n * sizeof(double));
    if (L == NULL)
    {
        fprintf(stderr, "Memory allocation for array L failed! \n");
        return 1;
    }

    // Allocate memory slots for the U-matrix.
    double *U = malloc(m * n * sizeof(double));
    if (U == NULL)
    {
        fprintf(stderr, "Memory allocation for array U failed! \n");
        return 1;
    }

    // Allocate memory slots for the y-vector.
    double *y = malloc(m * sizeof(double));
    if (y == NULL)
    {
        fprintf(stderr, "Memory allocation for array y failed! \n");
    }

    free(L);
    free(U);
    free(y);

    L = NULL;
    U = NULL;   
    y = NULL;

    return 0;
}

int init_L_matrix (double *A, int m, double *L, double *U)
{
    return 0;
}

int init_U_matrix (double *A, int m, double *L, double *U)
{
    return 0;
}


