#ifndef ELEMENTS_H
#define ELEMENTS_H

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

int initQuadElement (quadElement *element, int id, int gp, int nnodesElement, int dof, int *nodeids, double *coords);
void commitTrialValuesAtGaussPoints (quadElement *element, int nelements);
int displacementStrain (double *epsilon, double *B, double *ue, int nnodesElement, int dof);


#endif