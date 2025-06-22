#ifndef ELEMENTS_H
#define ELEMENTS_H

#include <stdio.h>
#include <stdlib.h>

#include "matrixArithmetic.h"

typedef struct { // REMOVE: This is moved to 'materials.c'
    int id;
    double E;
    double v;
    double H;
    double yield;
} material;


typedef struct { // REMOVE: This is moved to 'nodes.c'
    int id;
    double x;
    double y;
    double z;
    int material;
    int dof;
    int dim;
} node;


typedef struct { // REMOVE: This will be generalized to handle any element type and size.

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

typedef struct 
{
    // Element struct
    int id;
    int gp;
    int nnodes;
    int dof;
    int dim;
    int mat;
    int *nodeids;

    // Accumulated stress and strain
    double *sigma;
    double *epsilonP;
    double *epsilonBarP;

    // Trial stress and strain
    double trialSigma;
    double trialEpsilonP;
    double trialEpsilonBarP;
} element;

void commitTrialValuesAtGaussPoints (quadElement *element, int nelements); // REMOVE once the new element struct and helper functions are used.
int elementInit (element e, int gp, int dof, int dim, int nnodes, int mat, int flag);
int elementSetElasticTrialValues (element e, int gp, double sigma);
int elementSetPlasticTrialValues (element e, int gp, double epsilonP, double epsilonBarP);
void elementCommitElasticTrialValues (element e);
void elementCommitPlasticTrialValues (element e);


#endif