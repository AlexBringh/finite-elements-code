#ifndef FORCE_VECTOR_H
#define FORCE_VECTOR_H

void initElementInternalForceVector (double *FintE, int Kem);
void elementInternalForceVector (double *FintE, double *Btrans, double *sigma, double detJ, double weight, int Kem, int Bn);
void initGlobalInternalForceVector (double *Fint, int Km);
void globalInternalForceVector (double *Fint, double *FintE, int *nodeids, int dof, int nnodesElement, int Km, int Kem);

#endif