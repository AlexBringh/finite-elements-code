#ifndef SKYLINE_H
#define SKYLINE_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define SKYLINE_BUFFER 16  // Additional space when reallocating

typedef struct
{
    int n;              // Matrix dimension
    double *cellData;   // Non-zero values (skyline storage)
    int *colIndex;      // Index into cellData where each column starts
    int *colTop;        // Row index of the top non-zero entry in each column
    int storedCells;    // Actual number of stored values in cellData
    int capacity;       // Capacity allocated to cellData
} skylineMatrix;

void initSkylineMatrix(skylineMatrix *mat, int n);
void resetSkylineMatrix (skylineMatrix *mat);
void addToSkyline(skylineMatrix *mat, int i, int j, double value);
double getFromSkyline(skylineMatrix *mat, int i, int j);
void solveSkylineSystem(skylineMatrix *mat, double *r, double *u);

#endif