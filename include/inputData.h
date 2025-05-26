#ifndef INPUT_DATA_H
#define INPUT_DATA_H

#include "elements.h"

int loadMeshSize (char *path, int *nelements, int *nnodes, int *nnodesElement);
int loadMaterialData (char *path, double *E, double *v, double *sigmaYieldInitial, double *H);
int loadGeometry (char *path, quadElement *element, int *nelements);
int inputData (quadElement *element, int nElements, int nnodesElement);

#endif  