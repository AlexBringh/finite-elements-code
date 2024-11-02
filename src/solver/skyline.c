//
//	Author: Alexander B. Ringheim
//	Project: Finite Element Method solver for Plasticity analysis of materials and structures.
//	Date of creation: 01.11.2024
//

#include "skyline.h"


skylineMatrix* initSkylineMatrix (int n)
{
    // TODO: Documentation

    skylineMatrix *matrix = malloc(sizeof(skylineMatrix));

    matrix->n = n;
    matrix->cellData = malloc((n * (n + 1) / 2) * sizeof(double));
    if (matrix->cellData == NULL)
    {
        printf("Error (skyline): Memory allocation for array 'cellData' failed! \n");
        exit(1);
    }
    matrix->colIndex = malloc((n + 1) * sizeof(int));
    if(matrix->colIndex == NULL)
    {
        printf("Error (skyline): Memory allocation for array 'colIndex' failed! \n");
        exit(1);
    }

    for (int i = 0; i < (n + 1); i++)
    {
        *(matrix->colIndex + i) = 0;
    }

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

int updateSkylineElement (skylineMatrix* matrix, int m, int n, double val)
{
    // TODO: Documentation

    // TODO: Consider if this function is necessary and if so, implement it. If not, remove it.
    return 0;
}
