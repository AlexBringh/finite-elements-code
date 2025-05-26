#ifndef STIFFNESS_MATRIX_H
#define STIFFNESS_MATRIX_H

#include "elements.h"
#include "skyline.h"

int quadBMatrix (double *B, double *Btrans, double *Jinv, double *NiPxi, double *NiPeta);
int elementStiffnessMatrix (double *Ke, int nnodesElement, int DOF, double *B, double *Btrans, double *D, double detJ, double weight);
int initElementStiffnessMatrix (double *Ke, int Kem);
int elasticDMatrixPlaneStress (double *D, double E, double v);
int elasticDMatrixPlaneStrain (double *D, double E, double v);
int initGlobalStiffnessMatrix (double *K, int Km);
int globalStiffnessMatrix (double *K, double *Ke, int *nodeids, int dof, int nnodesElement, int Km, int Kem);
void applyFixedDisplacementStiffnessMatrix (int *uFixed, double *K, int Km);
int globalSkylineStiffnessMatrix (skylineMatrix *K, double *Ke, int *nodeids, int dof, int nnodesElement, int Kem);
void applyFixedDisplacementSkylineStiffnessMatrix (int *uFixed, skylineMatrix *K);

#endif