#ifndef POST_PROCESSING_H
#define POST_PROCESSING_H

#include "elements.h"
#include "shapeFunctions.h"

void postProcessingElastoPlastic (double *u, quadElement *elements, shapeFunctions2D *sf, int m);

#endif