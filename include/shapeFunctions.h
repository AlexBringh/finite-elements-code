#ifndef SHAPE_FUNCTIONS_H
#define SHAPE_FUNCTIONS_H

int quadShapeFunction (double *Ni, double *NiPxi, double *NiPeta);
int triangularShapeFunction ();
int hexShapeFunction ();
int tetraShapeFunction ();

#endif