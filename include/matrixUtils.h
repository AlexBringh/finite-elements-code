#ifndef MATRIX_UTILS_H
#define MATRIX_UTILS_H

#include <stdio.h>
#include <stdlib.h>

void printMatrix (double *A, int m, int n);
void printPreciseMatrix (double *A, int m, int n);
int augmentedMatrix (double *Ab, double *A, double *B, int m);

#endif