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

    // Store the values for the column index.
    for (int i = 0; i < (n + 1); i++)
    {
        *(matrix->colIndex + i) = *(startRow + i);
    }

    /* Commented out because this is likely redundant. // TODO: Remove
    // Transfer the non-zero elements from the A-matrix to the skyline matrix.
    int c = 0; // Counter for address positioning of cellData
    for (int i = 0; i < n + 1; i++) // Loop through the columns
    {
        for (int j = *(matrix->colIndex + i); j <= i; j++) // Start at the stored starting index for the column, and loop through the rows till the diagonal is reached.
        {
            *(matrix->cellData + c) = *(A + j*m + i); // The column values are stored next to each other in memory, rather than the row values being next to each other.
            c++;
        }
    }*/

    return matrix;
}

int addSkylineElement (skylineMatrix* matrix, int m, int n, double val)
{
    // TODO: Documentation

    int start = (matrix->colIndex + n);
    int stop = (matrix->colIndex + (n + 1));
    int pos = start + m;

    if (start <= pos && pos < stop)
    {
        *(matrix->cellData + pos) = val;
        return 0;
    }
    else 
    {
        printf("Error (skyline): Position of skyline profile out of range! \n");
        return 1;
    }
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
