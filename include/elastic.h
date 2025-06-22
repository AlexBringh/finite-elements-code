#ifndef ELASTIC_H
#define ELASTIC_H

#include <stdio.h>
#include <stdlib.h>

#include "elements.h"
#include "nodes.h"
#include "materials.h"

int displacementStrain (double *epsilon, double *B, double *ue, int nnodesElement, int dof);

#endif