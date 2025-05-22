//
//	Author: Alexander B. Ringheim
//	Project: Finite Element Method solver for Plasticity analysis of materials and structures.
//	Date of creation: 01.11.2024
//

#include "skyline.h"
#include <stdio.h>
#include <stdlib.h>

// Function prototypes only used for internal functionality that should be hidden from others.
//int findSkylineStartingRow (int *startRow, double *A, int m);


skylineMatrix* initSkylineMatrix (int n, int *startRow)
{
    // TODO: Documentation

    // Assuming that the matrix will always be cubic, number of rows, m, is set equal to number of columns, n.
    int m = n;

    // Allocate memory for the skyline matrix structure.
    skylineMatrix *matrix = malloc(sizeof(skylineMatrix));
    if (matrix == NULL)
    {
        fprintf(stderr, "Error (skyline): Memory allocation for skyline matrix failed! \n");
        exit(1);
    }

    // Store number of columns in the skyline matrix structure.
    matrix->n = n;

    // Calculate the number of elements in the skyline matrix. Assuming that the matrix is symmetric about the diagonal. Adds 50% buffer for fill-ins during crout reduction later.
    int elements = (m + 1) * n / 2 * 1.5;

    // Allocate memory for the cell data of the skyline matrix.
    matrix->cellData = malloc(elements * sizeof(double));
    if (matrix->cellData == NULL)
    {
        fprintf(stderr, "Error (skyline): Memory allocation for array 'cellData' failed! \n");
        exit(1);
    }

    // Allocate memory for the column index.
    matrix->colIndex = malloc((n + 1) * sizeof(int));
    if(matrix->colIndex == NULL)
    {
        fprintf(stderr, "Error (skyline): Memory allocation for array 'colIndex' failed! \n");
        exit(1);
    }

    matrix->colTop = malloc(n * sizeof(int));
    if (matrix->colTop == NULL)
    {
        fprintf(stderr, "Error (skyline): Memory allocation for array 'colTop' failed! \n");
        exit(1);
    }

    // Store the values for the column index.
    for (int i = 0; i < (n + 1); i++)
    {
        *(matrix->colTop + i) = *(startRow + i);
    }

    *(matrix->colIndex) = 0;
    for (int i = 0; i < n; i++)
    {
        *(matrix->colIndex + i + i) = *(matrix->colIndex + i) + (i - *(startRow + i) + 1);
    }

    return matrix;
}

int addSkylineElement (skylineMatrix* matrix, int m, int n, double val)
{
    /*
    Adds data values to the cells of the skyline matrix, based upon the starting values of the colTop values.
    */

   // Check that the data is within the skyline structure bounds.
    if (m < *(matrix->colTop + n) || m > n)
    {
        fprintf(stderr, "Error (skyline): Could not add element to skyline matrix. Index out of range!");
        return 1;
    }

    // Sett data at index.
    int index = *(matrix->colIndex + n) + (m - *(matrix->colTop + n));
    *(matrix->cellData + index) = val; 
}

int setSkylineElement (skylineMatrix* matrix, int m, int n, double val)
{
    // TODO: Documentation

    // Check that the index is valid
    if (m > n) 
    {
        fprintf(stderr, "Error (skyline): Could not set element in skyline matrix. Index out of range!");
        return 1;
    }

    // Calculate the target index in the 'data' array.
    int colOffset = m - *(matrix->colTop + n);
    int index = *(matrix->colIndex + n) + colOffset;

    // Check that the index pos is within the skyline matrix structure
    if (m < *(matrix->colTop + n))
    {
        // If the row, m, is outside of the skyline rows for this column, shift 'colTop' and reassign the 'colIndex' values.
        int shift = *(matrix->colTop + n) - m;
        // Shift 'colTop'
        *(matrix->colTop + n) = m; 
        // Reassign / add to colIndex values.
        for (int k = n; k <= matrix->n; k++)
        {
            *(matrix->colIndex + k) += shift;   
        }

        index = *(matrix->colIndex + n);

        // Reallocate and shift 'data' array.
        matrix->cellData = realloc(matrix->cellData, sizeof(double) * *(matrix->colIndex + (matrix->n)) - *(matrix->colIndex));

        // Shift the elements in the 'data' array.
        for (int k = *(matrix->colIndex + (matrix->n)); k > index; k--)
        {
            *(matrix->cellData + k) = *(matrix->cellData + k - shift);
            if (k < index + shift)
            {
                *(matrix->cellData + k - shift) = 0; // Init data cell with 0.
            }
        }
    }

    // Insert / update the value for the given cell.
    *(matrix->cellData + index) = val;
}

double getSkylineElement (skylineMatrix* matrix, int m, int n)
{
    // TODO: Documentation

    int start = *(matrix->colIndex + n);
    int stop = *(matrix->colIndex + (n + 1));
    int pos = start + m;

    if (start <= pos && pos < stop)
    {
        return *(matrix->cellData + pos);
    }
    else
    {
        return 0;
    }
}

/* Commented out because this is likely redundant. // TODO: Remove
int findSkylineStartingRow (int *startRow, double *A, int m)
{
    // TODO: Documentation
    int n = m;

    for (int i = 0; i < n; i++) // Columns to check through
    {
        int c = 0; // Counter
        for (int j = 0; j < m; j++) // Check through the rows, end once a non-zero member is found.
        {
            if (*(A + j*n + i) == 0) // If the element at the j-th row on the i-th column is 0, add to the counter.
            {
                c++;
            }
            else // If the element at the j-th row on the i-th column is NOT zero, then we have our starting-point for the column and the loop can be ended.
            {
                break; 
            }
        }

        *(startRow + i) = c; // Add the counter variable for the current column and continue to the next column.
    }
}*/

int resetSkylineMatrix (skylineMatrix *matrix)
{
    
}
