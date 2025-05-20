#ifndef ELEMENTS_H
#define ELEMENTS_H

typedef struct {

    int id; // Element id
    int gp; // Number of Gauss Points
    int nnodes; // Number of nodes
    int dof; // Degrees of freedom per node
    int *nodeids; // Node connectivity numbering
    int *coords; // Global coordinates of connected nodes.

    // Accumulated / commited stress and plastic strain.
    double *sigma;
    double *epsilonP;
    double epsilonBarP;

    // Trial stress and strain
    double *trialSigma;
    double *trialEpsilonP;
    double trialEpsilonBarP;

} quadElement;

int displacementStrain (double *epsilon, double *B, double *ue, int nnodesElement, int dof);


#endif