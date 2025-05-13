#ifndef STIFFNESS_MATRIX_H
#define STIFFNESS_MATRIX_H

int quadBMatrix (double *B, double *Btrans, double *Jinv, double *NiPxi, double *NiPeta);
int quadElementStiffnessMatrix (double *Ke, int gp, int DOF, double *B, double *Btrans, double *D, double detJ, double *weights);

#endif