#ifndef ELEMENTS_H
#define ELEMENTS_H

typedef struct {

    int id; // Element id
    int nnodes; // Number of nodes
    int dof; // Degrees of freedom per node
    int *connectivity; // node connectivity numbering
    int *coords; // Coordinates of connected nodes from local numbering 1 to nnodes.
    int gp; // Number of Gauss Points

    double epsilon_p;
    double epsilon_bar_p;
    double sigma;

} quadElement;

#endif