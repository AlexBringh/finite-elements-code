//
//	Author: Alexander B. Ringheim
//	Project: Finite Element Method solver for Plasticity analysis of materials and structures.
//	Date of creation: 01.11.2024
//

#include "skyline.h"
#include <stdio.h>
#include <stdlib.h>

// Local function prototypes for functions that are only to be used inside this file.
double* allocateAndZero(int size);
void ensureCapacity(skylineMatrix *mat, int extraNeeded);

// Initialize the skyline matrix
double* allocateAndZero(int size) 
{
    double* ptr = (double*)calloc(size, sizeof(double));
    if (!ptr) 
    {
        fprintf(stderr, "Failed to allocate memory.\n");
        exit(EXIT_FAILURE);
    }
    return ptr;
}

void initSkylineMatrix(skylineMatrix *mat, int n) 
{
    mat->n = n;
    mat->colTop = (int*)calloc(n, sizeof(int));
    mat->colIndex = (int*)calloc(n + 1, sizeof(int)); // n+1 to simplify last column
    mat->capacity = n * SKYLINE_BUFFER;
    mat->cellData = allocateAndZero(mat->capacity);
    mat->storedCells = 0;
    for (int i = 0; i <= n; i++) {
        mat->colIndex[i] = 0;
    }
}

// Reset all stored values to zero
void resetSkylineMatrix(skylineMatrix *mat) 
{
    memset(mat->cellData, 0, sizeof(double) * mat->storedCells);
}

// Expand storage if needed
void ensureCapacity(skylineMatrix *mat, int extraNeeded) 
{
    if (mat->storedCells + extraNeeded > mat->capacity) 
    {
        mat->capacity += extraNeeded + SKYLINE_BUFFER;
        mat->cellData = (double*)realloc(mat->cellData, sizeof(double) * mat->capacity);
        if (!mat->cellData) 
        {
            fprintf(stderr, "Failed to reallocate memory.\n");
            exit(EXIT_FAILURE);
        }
    }
}

// Add a value to K[i][j] += value
void addToSkyline(skylineMatrix *mat, int i, int j, double value) {
    if (i > j) return;  // only store upper triangle

    int col = j;
    if (i < mat->colTop[col]) {
        int move = mat->colTop[col] - i;

        // Compute offset to start of column 'col'
        int startOffset = 0;
        for (int k = 0; k < col; ++k) {
            startOffset += (k - mat->colTop[k] + 1);
        }

        int colHeight = col - mat->colTop[col] + 1;
        int colLength = colHeight;

        ensureCapacity(mat, move);

        // Shift memory to make room
        memmove(mat->cellData + startOffset + move, mat->cellData + startOffset, sizeof(double) * colLength);
        memset(mat->cellData + startOffset, 0, sizeof(double) * move);

        mat->colTop[col] = i;
        mat->storedCells += move;

        for (int k = col + 1; k <= mat->n; ++k) {
            mat->colIndex[k] += move;
        }
    }

    // Find offset again now that top may have changed
    int startOffset = 0;
    for (int k = 0; k < j; ++k) {
        startOffset += (k - mat->colTop[k] + 1);
    }
    int offset = j - i;
    int idx = startOffset + offset;
    mat->cellData[idx] += value;
}  

// Retrieve a value from the matrix
double getFromSkyline(skylineMatrix *mat, int i, int j) {
    if (i > j) {
        int temp = i; i = j; j = temp;
    }

    if (i < mat->colTop[j]) {
        return 0.0;
    }

    int startOffset = 0;
    for (int k = 0; k < j; ++k) {
        startOffset += (k - mat->colTop[k] + 1);
    }

    int offset = j - i;
    int idx = startOffset + offset;
    return mat->cellData[idx];
}

// Crout reduction solver (no pivoting)
void solveSkylineSystem(skylineMatrix *mat, double *r, double *u) 
{
    int n = mat->n;
    double *y = (double*)calloc(n, sizeof(double));

    // Forward substitution (Ly = r)
    for (int i = 0; i < n; i++) 
    {
        double sum = 0.0;
        for (int j = mat->colTop[i]; j < i; j++) {
            sum += getFromSkyline(mat, i, j) * y[j];
        }
        y[i] = (r[i] - sum) / getFromSkyline(mat, i, i);
    }

    // Backward substitution (Lu = y)
    for (int i = n - 1; i >= 0; i--) 
    {
        double sum = 0.0;
        for (int j = i + 1; j < n; j++) 
        {
            if (i >= mat->colTop[j]) 
            {
                sum += getFromSkyline(mat, j, i) * u[j];
            }
        }
        u[i] = (y[i] - sum) / getFromSkyline(mat, i, i);
    }

    free(y);
}

