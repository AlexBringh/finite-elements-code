#ifndef INPUT_DATA_H
#define INPUT_DATA_H

#include "elements.h"

int readCSV(const char *filename, node **nodes_out, int *nnodes_out, int *nnodesElement, int *nelements, int *gp, material **materials_out, int *nmaterials_out, int **uFixed_out, double **fApplied_out);
int readElements(const char *filename, quadElement **elements_out, int *nelements_out, node *nodes, int nnodes, int gp, int nnodes_per_elem, int dof);

#endif  