//
//	Author: Alexander B. Ringheim
//	Project: Finite Element Method solver for Plasticity analysis of materials and structures.
//	Date of creation: 15.11.2024
//

#include "solveCroutSkyline.h"
#include "matrixUtils.h"
#include "printUtils.h"
#include <stdio.h>
#include <stdlib.h>

// Local function prototypes
int luDecompositionSkyline (skylineMatrix *matrix, double *L, double *U);

int croutSkyline (skylineMatrix *matrix)
{
    // TODO: DOCUMENTATION

    int n = matrix->n;
    int d = 30; // Constant used for printing the right amount of dashed lines.

    // Allocate memory slots for the L-matrix.
    double *L = malloc(n * n * sizeof(double));
    if (L == NULL)
    {
        fprintf(stderr, "Error (solveCroutSkyline): Memory allocation for array L failed! \n");
        printDashedLines(d);
        return 1;
    }

    // Allocate memory slots for the U-matrix.
    double *U = malloc(n * n * sizeof(double));
    if (U == NULL)
    {
        fprintf(stderr, "Error (solveCroutSkyline): Memory allocation for array U failed! \n");
        printDashedLines(d);
        return 1;
    }

    // Run LU decomposition
    if (luDecompositionSkyline(matrix, L, U))
    {
        fprintf(stderr, "Error (solveCroutSkyline): Could not perform LU-decomposition! \n");
        return 1;
    }
}

int luDecompositionSkyline (skylineMatrix *matrix, double *L, double *U)
{
    // TODO: DOCUMENTATION

    int n = matrix->n;

    for (int i = 0; i < n; i++)
    {
        int index = *(matrix->colIndex + i) + (i - *(matrix->colTop + i));
        *(U + index) = 1; // Diagonal of U is always 1.

        // Compute values for L on the i-th row.
        for (int j = *(matrix->colTop + i); j <= i; j++)
        {
            double sum = 0;
            for (int k = *(matrix->colTop + i); k < j; k++)
            {
                sum += getSkylineElement(matrix, j, k) * getSkylineElement(matrix, k, i);
            }
            int indexL = *(matrix->colIndex + i) + (j - *(matrix->colTop + i));
            *(L + indexL) = getSkylineElement(matrix, j, i) - sum;
        }

        for (int j = i + 1; j <= n; j++)
        {
            if (i < *(matrix->colTop + j))
            {
                continue;
            }

            double sum = 0;
            for (int k = *(matrix->colTop + j); k < i; k++)
            {
                sum += getSkylineElement(matrix, i, k) * getSkylineElement(matrix, k, j);
            }
            int indexU = *(matrix->colIndex + i) + (i - *(matrix->colTop + j));
            *(U + indexU) = (getSkylineElement(matrix, i, j) - sum) / getSkylineElement(matrix, i, i);
        }
    }
}