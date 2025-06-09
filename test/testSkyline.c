#include "testSkyline.h"
#include "skyline.h"
#include <stdio.h>

void testSkyline();

int main() 
{
    printf("Test skyline matrix: \n");
    skylineMatrix mat;
    int n = 3;
    initSkylineMatrix(&mat, n);

    // Manually add a small symmetric matrix:
    // [ 4  -1   0 ]
    // [-1   4  -1 ]
    // [ 0  -1   4 ]
    addToSkyline(&mat, 0.0, 0.0, 4.0);
    addToSkyline(&mat, 0.0, 1.0, -2.0);
    addToSkyline(&mat, 1.0, 1.0, 5.0);
    addToSkyline(&mat, 1.0, 2.0, -1.0);
    addToSkyline(&mat, 2.0, 2.0, 8.0);
    addToSkyline(&mat, 0.0, 2.0, 3.0);
    addToSkyline(&mat, 0.0, 1.0, 7.0);
    addToSkyline(&mat, 2.0, 2.0, 2.0);

    printf("Matrix: \n");
    for(int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            printf("%.4f   ", getFromSkyline(&mat, i, j));
        }
        printf("\n");
    }

    double r[3] = {1.0, 2.0, 3.0};
    double u[3] = {0.0, 0.0, 0.0};

    solveSkylineSystem(&mat, r, u);

    printf("Displacement vector u:\n");
    for (int i = 0; i < n; i++) {
        printf("u[%d] = %f\n", i, u[i]);
    }

    printf("Stored cells: %d\n", mat.storedCells);

    resetSkylineMatrix(&mat);

    printf("Reset Matrix: \n");
    for(int i = 0; i < mat.n; i++)
    {
        for (int j = 0; j < mat.n; j++)
        {
            printf("%.4f   ", getFromSkyline(&mat, i, j));
        }
        printf("\n");
    }

    free(mat.cellData);
    free(mat.colIndex);
    free(mat.colTop);
    return 0;
}


void testSkyline()
{
    // TODO testSkyline
}