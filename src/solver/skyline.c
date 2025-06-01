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

double* allocateAndZero(int size) 
{
    /*
        Helper function for allocating memory

        Inputs:
        int size    -> Size of memroy chunk to be allocated

        Output: 
        double *ptr -> Pointer with allocated memory.
    */

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
    /*
        Initializes and allocates memory to the matrix, and sets it's size, n.

        Inputs:
        skylineMatrix *mat  -> Pointer to skyline matrix.
        int n               -> Size, n, of the matrix.
    */

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


void resetSkylineMatrix(skylineMatrix *mat) 
{
    /*
        Resets the skyline matrix by setting all the stored cell data to 0.

        Input:
        skylineMatrix *mat  -> Pointer to skyline matrix.
    */
    memset(mat->cellData, 0, sizeof(double) * mat->storedCells);
}


void ensureCapacity(skylineMatrix *mat, int extraNeeded) 
{
    /*
        Ensures there is enough capacity in the matrix' cell storage.
        If not, allocates more memory.
        Local function, not to be used outside of this file.

        Inputs:
        skylineMatrix *mat  -> Pointer to skyline matrix.
        int extraNeeded     -> How many memory chunks are needed.
    */

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


void addToSkyline(skylineMatrix *mat, int i, int j, double value) 
{
    /*
        Adds new value to the skyline matrix or adds to existing ones.
        Checks first row-index i > column-index j. 
        If so does nothing, and assumes the symmetric index will come 
        later or has come earlier.
        Checks if the value sent in is 0, if so does not need to store it.
        Checks if there is an existing value stored for the indices.
        If so, adds the value to that.
        If not, the structure must be updated to make room for the new value.
        Calculates how far to move the memory over, checks if there is enough
        allocated memory.
        If not, allocates more memroy, alongside an extra buffer for later,
        and moves the existing memory over. Sets the freed slots to 0,
        then adds the value to the correct cell. Updates the column top indexes.

        Inputs:
        skylineMatrix *mat  -> Pointer to the skyline matrix.
        int i               -> Row-index, i.
        int j               -> Column index, j.
        double value        -> Value to be stored int he matrix.
    */

    if (i > j) return;  // only store upper triangle
    if (value == 0) return; // Never store 0, zeroes are implied by the lack of nonzero values for the skyline matrix method.

    int col = j;
    if (i < mat->colTop[col]) 
    {
        int move = mat->colTop[col] - i;

        // Compute offset to start of column 'col'
        int startOffset = 0;
        for (int k = 0; k < col; ++k) 
        {
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

        for (int k = col + 1; k <= mat->n; ++k) 
        {
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


double getFromSkyline(skylineMatrix *mat, int i, int j) 
{
    /*
        Gets cells at row-index i and column-index j.
        If the value is not stored, returns 0.
        If the row index i > j, then swap the indices and take advantage of symmetry.

        Inputs: 
        skylineMatrix *mat  -> Pointer to skyline matrix struct.
        int i               -> Row index, i
        int j               -> Column index, j

        Output:
        double value        -> Returns the value at the cell or 0 if not stored.
    */

    // Use symmetry, if the row index is greater than the column index, swap them.
    if (i > j) 
    {
        int temp = i; i = j; j = temp;
    }

    // If the row index is larger than the top column index, return 0.
    if (i < mat->colTop[j]) 
    {
        return 0.0;
    }

    // Find the offset from the start to the column, j.
    int startOffset = 0;
    for (int k = 0; k < j; ++k) 
    {
        startOffset += (k - mat->colTop[k] + 1);
    }

    int offset = j - i;
    int idx = startOffset + offset;
    return mat->cellData[idx];
}


void solveSkylineSystem(skylineMatrix *mat, double *r, double *u) 
{
    /*
        Solves the skyline matrix using Crout reduction without pivoting.
        Good luck to anyone trying to implement pivoting with the skyline method.
        To avoid large memory usage, does not make dedicated L and U matrices,
        does the calculations directly on the y and u matrices (u being x in Lx = y).

        Input:
        skylineMatrix *mat  -> Pointer to skyline matrix struct
        double *r           -> Pointer to residual force vector
        double *u           -> Pointer to displacement vector, u. Results are stored here.
    */

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
