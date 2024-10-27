//
//	Author: Alexander B. Ringheim
//	Project: Finite Element Method solver for Plasticity analysis of materials and structures.
//	Date of creation: 21.10.2024
//

#include "croutReduction.h"
#include "matrixUtils.h"

// Function prototypes that are limited access to only this file.
void partialPivot (double *A, int m, double *P, int k);
int luDecomposition (double *A, int m, double *L, double *U, double *P);
int permutateVectorB (double *B, double *P, int m);
int forwardSubstitution (double *B, int m, double *L, double *U, double *y);
int backwardSubstitution (double *x, double *y, double *U, int m);

int croutReduction (double *A, int m, double *x, double *B)
{
    // TODO documentation
    int n = m; // n represents the columns of the A-matrix.

    // Allocate memory slots for the L-matrix.
    double *L = malloc(m * n * sizeof(double));
    if (L == NULL)
    {
        fprintf(stderr, "Error (croutReduction): Memory allocation for array L failed! \n");
        return 1;
    }

    // Allocate memory slots for the U-matrix.
    double *U = malloc(m * n * sizeof(double));
    if (U == NULL)
    {
        fprintf(stderr, "Error (croutReduction): Memory allocation for array U failed! \n");
        return 1;
    }

    // Allocate memory slots for the y-vector.
    double *y = malloc(m * sizeof(double));
    if (y == NULL)
    {
        fprintf(stderr, "Error (croutReduction): Memory allocation for array y failed! \n");
        return 1;
    }

    // Allocate memory slots for the permutation (P) vector.
    double *P = malloc(m * sizeof(double));
    if (P == NULL)
    {
        fprintf(stderr, "Error (croutReduction): Memory allocation for array P failed! \n");
        return 1;
    }
    // Init the permutation vector
    for (int i = 0; i < m; i++)
    {
        *(P + i) = i;
    }

    // Print out the input A-matrix and B-vector
    printf("---------------------------------------------\n");
    printf("Crout redcution: \n\nInputed A-matrix: \n");
    printMatrix(A, m, n);
    printf("Inputed B-matrix: \n");
    printMatrix(B, m, 1);

    luDecomposition(A, m, L, U, P);
    permutateVectorB(B, P, m);

    printf("Result of the system: \n");
    //TODO: UNCOMMENT: printMatrix(x, m, 1);
    printf("---------------------------------------------\n\n");

    free(L);
    free(U);
    free(y);
    free(P);

    L = NULL;
    U = NULL;   
    y = NULL;
    P = NULL;

    return 0;
}

void partialPivot (double *A, int m, double *P, int k)
{
    /*
        Partial pivots the A-matrix and correspondingly the B-vector so that the largest values for the i-th column will go down diagonally. 
        If the larges value of x2 belongs to the row that also has the largest value of x1, then the value for x1 will be favored so that x1 gets its largest value on the diagonal. 
        x2 can then have its largest value, excluding the value on the same row as x1, on the next row on the diagonal.

        This process does, unfortuantely add complexity and therefore computation time, however it is absolutely necessary to limit the amount of computational error or algorithm breakdowns caused by potential small values being placed on the diagonal of the L-matrix later.
    */

    int n = m;

    int maxIndex = k;
    double maxVal = fabs(*(A + k * n + k));
    for (int i = k + 1; i < n; i++) {
        double val = fabs(*(A + i * n + k));
        if (val > maxVal) {
            maxIndex = i;
            maxVal = val;
        }
    }
    
    if (maxIndex != k) {
        for (int j = 0; j < n; j++) {
            double temp = *(A + k * n + j);
            *(A + k * n + j) = *(A + maxIndex * n + j);
            *(A + maxIndex * n + j) = temp;
        }
        int temp = *(P + k);
        *(P + k) = *(P + maxIndex);
        *(P + maxIndex) = temp;
    }
}

int luDecomposition (double *A, int m, double *L, double *U, double *P)
{
    // TODO documentation
    int n = m; // n represents columns.

    // Initialize L and U as matrices of 0s.
    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < n; j++)
        {
            *(L + i*n + j) = 0;
            *(U + i*n + j) = 0;
        }
        *(U + i*n + i) = 1;
    }

    for (int i = 0; i < m; i++) // i represents the row multiplier for the new matrices
    {
        // Partial pivoting at each iteration of the LU-decomposition.
        partialPivot(A, m, P, i);

        // Compute values for L on the i-th row.
        for (int j = i; j < n; j++) // j represents the column multiplier for the new matrices.
        {
            double lu = 0;

            for (int k = 0; k < i; k++)
            {
                lu += *(L + j*n + k) * *(U + k*n + i);
            }
            *(L + j*n + i) = *(A + j*n + i) - lu;
        }

        // Compute values for U on the i-th row.
        for (int j = i; j < n; j++)
        {
            if (*(L + i*n + i) == 0) // In case diagonal value of L_ii is 0, exit function and print an error.
            {
                fprintf(stderr, "Error (croutReduction): Could not complete LU-decomposition! Diagonal value of L is 0. \n");
                return 1;
            }

            double lu = 0;

            for (int k = 0; k < i; k++)
            {
                lu += *(L + i*n + k) * *(U + k*n + j);
            }
            *(U + i*n + j) = (*(A + i*n + j) - lu) / (*(L + i*n + i));
        }
        
        // The diagonal on the U-matrix is always 1.
        *(U + i*n + i) = 1; 
    }

    // Print the L- and U- matrices, primarily for debugging purposes.
    printf("L-matrix: \n");
    printMatrix(L, m, n);
    printf("U-matrix: \n");
    printMatrix(U, m, n);
    printf("\n");

    return 0;
}

int permutateVectorB (double *B, double *P, int m)
{
    // TODO documentation
    // Reserve memory for a temporary array for storing sorted B-values.
    double *tempB = malloc(m * sizeof(double));
    if (tempB == NULL)
    {
        fprintf(stderr, "Error (croutReduction): Memory allocation for array tempB failed! \n");
        return 1;
    }

    // Store values from B into correct slots in tempB.
    for (int i = 0; i < m; i++)
    {
        int iterator = *(P + i); // P stores values from 0 to (m - 1), corresponding to the space that the row was moved to during partial pivoting. This iterator shows what value from B the i-th row should have.
        *(tempB + i) = *(B + iterator);
    }

    // Store the values from tempB into the B-array.
    for (int i = 0; i < m; i++)
    {
        *(B + i) = *(tempB + i);
    }

    // Free the memory and remove the pointer to tempB.
    free(tempB);
    tempB = NULL;

    return 0;
}

int forwardSubstitution (double *B, int m, double *L, double *U, double *y)
{
    // TODO documentation
}

int backwardSubstitution (double *x, double *y, double *U, int m)
{
    // TODO documentation
}

