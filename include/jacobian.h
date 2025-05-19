#ifndef JACOBIAN_H
#define JACOBIAN_H

int quadJacobian (double *J, double *Jinv, int m, double *detJ, double *Ni, double *NiPxi, double *NiPeta, double *coords, int dof);

#endif