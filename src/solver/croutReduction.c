//
//	Author: Alexander B. Ringheim
//	Project: Finite Element Method solver for Plasticity analysis of materials and structures.
//	Date of creation: 21.10.2024
//

#include "croutReduction.h"
#include "matrixUtils.h"
#include "printUtils.h"

// Function prototypes that are limited access to only this file.
void partialPivot (double *A, int m, int *P, int k);
int luDecomposition (double *A, int m, double *L, double *U, int *P);
int permutateVectorB (double *B, int *P, int m);
int forwardSubstitution (double *B, int m, double *L, double *y);
int backwardSubstitution (double *x, double *y, double *U, int m);

int croutReduction (double *A, int m, double *x, double *B)
{
    /*
        Performs Crout Reduction with partial pivoting.

        Inputs:
        double *A -> Pointer to matrix holding the system of equations to be solved, of size m x m.
        double *x -> Pointer to vector holding the unknown values to be solved for, of size m x 1.
        double *B -> Pointer to vector holding the 
    */

    int n = m; // n represents the columns of the A-matrix.
    int d = 30; // Constant used for printing the right amount of dashed lines.

    // Allocate memory slots for the L-matrix.
    double *L = malloc(m * n * sizeof(double));
    if (L == NULL)
    {
        fprintf(stderr, "Error (croutReduction): Memory allocation for array L failed! \n");
        printDashedLines(d);
        return 1;
    }

    // Allocate memory slots for the U-matrix.
    double *U = malloc(m * n * sizeof(double));
    if (U == NULL)
    {
        fprintf(stderr, "Error (croutReduction): Memory allocation for array U failed! \n");
        printDashedLines(d);
        return 1;
    }

    // Allocate memory slots for the y-vector.
    double *y = malloc(m * sizeof(double));
    if (y == NULL)
    {
        fprintf(stderr, "Error (croutReduction): Memory allocation for array y failed! \n");
        printDashedLines(d);
        return 1;
    }

    // Allocate memory slots for the permutation (P) vector.
    int *P = malloc(m * sizeof(int));
    if (P == NULL)
    {
        fprintf(stderr, "Error (croutReduction): Memory allocation for array P failed! \n");
        printDashedLines(d);
        return 1;
    }
    // Init the permutation vector
    for (int i = 0; i < m; i++)
    {
        *(P + i) = i;
    }

    // Print out the input A-matrix and B-vector
    /*
    printDashedLines(d);
    printf("Crout redcution: \n\nInputed stiffness matrix: \n");
    printMatrix(A, m, n);
    printf("Inputed residual vector: \n");
    printMatrix(B, 1, m);
    */

    // Step 1: Decompose the A-matrix into L- and U-matrices (lower- and upper triangular matrices).
    luDecomposition(A, m, L, U, P);

    // Step 2: Permutate the B-vector to orient the vector in accordance with the pivoted L- and U- matrices.
    permutateVectorB(B, P, m);

    // Step 3: Perform forward substitution to obtain the y-vector values.
    forwardSubstitution(B, m, L, y);

    // Step 4: Perform backward substitution to obtain the x-vector values (the result of the system).
    backwardSubstitution(x, y, U, m);

    /*
    printf("Result of the system: \n");
    printPreciseMatrix(x, 1, m);
    printDashedLines(d);
    */

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

void partialPivot (double *A, int m, int *P, int k)
{
    /*
        Partial pivots the A-matrix and correspondingly the B-vector so that the largest values for the i-th column will go down diagonally. 
        If the larges value of x2 belongs to the row that also has the largest value of x1, then the value for x1 will be favored so that x1 gets its largest value on the diagonal. 
        x2 can then have its largest value, excluding the value on the same row as x1, on the next row on the diagonal.

        This process does, unfortuantely add complexity and therefore computation time, however it is absolutely necessary to limit the amount of computational error or algorithm breakdowns caused by potential small values being placed on the diagonal of the L-matrix later.

        Inputs:
        double *A   -> Pointer to matrix A. Results are stored here
        int m       -> Size m of matrix and vectors
        int *P      -> Partial pivoting indices
        int k       -> Column-index, k
    */

    // Check that the values in P are within bounds.
    for (int i = 0; i < m; i++)
    {
        if (*(P + i) > (m - 1) || *(P + i) < 0)
        {
            fprintf(stderr, "Error: P[%d] = %d out of bounds [0. %d]: ", i, *(P + i), (m - 1));
            exit(1);
        }
    }

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

int luDecomposition (double *A, int m, double *L, double *U, int *P)
{
    /*
        Decompose the matrix, A into upper- and lower diagonal matrices, U and L,
        such that LU = A. Performs partial pivoting along the steps to reduce
        the risk of division by 0 errors.

        Inputs:
        double *A   -> Pointer to matrix A, the system to be solved
        int m       -> Size m for the matrices.
        double *L   -> Pointer to matrix L
        double *U   -> Pointer to matrix U
        int *P      -> Pointer to pivot indices
    */
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
    /*
    printf("L-matrix: \n");
    printMatrix(L, m, n);
    printf("U-matrix: \n");
    printMatrix(U, m, n);
    printf("\n");
    */

    return 0;
}

int permutateVectorB (double *B, int *P, int m)
{
    /*
        Permutate vector B after partial pivoting to ensure the matrix indices match.
        
        Inputs:
        double *B   -> Pointer to B vector. Results are stored here.
        int *P      -> Pointer to pivot indices, P.
        int m       -> Size m of the vector, B.
    */
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

int forwardSubstitution (double *B, int m, double *L, double *y)
{
    /*
        Perform forward substitution Lb = y

        Inputs:
        double *B   -> Pointer to vector B
        int m       -> Size m of the matrix and vector
        double *L   -> Pointer to lower diagonal matrix
        double *y   -> Pointer to vector, y. Results are stored here.
    */
    int n = m;
    for (int i = 0; i < m; i++)
    {
        double ly = 0;
        for (int j = 0; j < i; j++)
        {
            ly += *(L + i*n + j) * *(y + j);
        }
        *(y + i) = ( *(B + i) - ly ) / *(L + i*n + i);
    }

    return 0;
}

int backwardSubstitution (double *x, double *y, double *U, int m)
{
    /*
        Perform backward substitution Lx = y to solve the system

        Inputs: 
        double *x   -> Pointer to the vector, x, where the solution is stored.
        double *y   -> Pointer to vector, y, from previous steps of the solution.
        double *U   -> Pointer to upper-diagonal matrix, U.
        int m       -> Size m of the matrix and vectors.
    */
    int n = m;
    for (int i = n - 1; i >= 0; i--)
    {
        double ux = 0;
        for (int j = i + 1; j < m; j++)
        {
            ux += *(U + i*n + j) * *(x + j);
        }
        *(x + i) = *(y + i) - ux;
    }

    return 0;
}
