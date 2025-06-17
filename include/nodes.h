#ifndef NODES_H
#define NODES_H

#include "materials.h"

typedef struct
{
    // Nodes struct
    int id;  // Node id (global)
    int dof; // Degrees of freedom
    int dim; // Dimensions
    int mat; // Material number
    double *coords;
} node;

#endif