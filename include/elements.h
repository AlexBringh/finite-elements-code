#ifndef ELEMENTS_H
#define ELEMENTS_H

#include <stdlib.h>

typedef struct {
    int id;
    double E;
    double v;
    double H;
    double yield;
} material;


typedef struct {
    int id;
    double x;
    double y;
    double z;
    int material;
    int dof;
    int dim;
} node;


typedef struct {

    int id; // Element id
    int gp; // Number of Gauss Points
    int nnodes; // Number of nodes
    int dof; // Degrees of freedom per node
    int *nodeids; // Node connectivity numbering
    double *coords; // Global coordinates of connected nodes.

    // Accumulated / commited stress and plastic strain.
    double *sigma; // Accumulated stress for all Gauss Points
    double *epsilonP; // Accumulated plastic strain tensors for all Gauss Points
    double *epsilonBarP; // Accumulated equivalent plastic strain for all Gauss Points

    // Trial stress and strain
    double *trialSigma; // 
    double *trialEpsilonP;
    double *trialEpsilonBarP;

} quadElement;

void commitTrialValuesAtGaussPoints (quadElement *element, int nelements);
int displacementStrain (double *epsilon, double *B, double *ue, int nnodesElement, int dof);


#endif