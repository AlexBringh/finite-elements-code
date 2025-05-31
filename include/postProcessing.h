#ifndef POST_PROCESSING_H
#define POST_PROCESSING_H

#include "elements.h"
#include "jacobian.h"
#include "shapeFunctions.h"
#include "stiffnessMatrix.h"
#include "vonMisesYield.h"


void postProcessingElastoPlastic (double *nodalStress, double *nodalStrain, double *nodalPlasticStrain, int nnodes, double *u, quadElement *elements, int nelements, shapeFunctions2D *sf, int m, double E);

#endif