#ifndef SKYLINE_H
#define SKYLINE_H

#include <stdio.h>
#include <stdlib.h>

typedef struct
{
    int n;
    double *cellData;
    int *colIndex;
} skylineMatrix;

skylineMatrix* initSkylineMatrix (int n, int *startRow);
int addSkylineElement (skylineMatrix* matrix, int m, int n, double val);
double getSkylineElement (skylineMatrix* matrix, int m, int n);

#endif