#ifndef STIFFNESS_MATRIX_H
#define STIFFNESS_MATRIX_H

int quadBMatrix (double *B, double *Btrans, double *Jinv, double *NiPxi, double *NiPeta);
int elementStiffnessMatrix (double *Ke, int nnodesElement, int DOF, double *B, double *Btrans, double *D, double detJ, double weight);
int initElementStiffnessMatrix (double *Ke, int nnodesElement, int DOF);
int elasticDMatrixPlaneStrain (double *D, double E, double v);

#endif