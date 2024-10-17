#ifndef GAUSS_ELIMINATION_H
#define GAUSS_ELIMINATION_H

#include <stdio.h>
#include <math.h>

int gauss_elimination (double *A, int n, double *B);
int forward_elimination (double *A, int n);
int backwards_substitution (double *A, int n, double *B);

#endif