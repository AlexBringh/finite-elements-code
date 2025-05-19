#ifndef ELEMENTS_H
#define ELEMENTS_H

typedef struct {

    int id; // Element id
    int gp; // Number of Gauss Points
    int nnodes; // Number of nodes
    int dof; // Degrees of freedom per node
    int *nodeids; // Node connectivity numbering
    int *coords; // Global coordinates of connected nodes.

    // Accumulated stress and plastic strain.
    double sigma;
    double epsilon_p;
    double epsilon_bar_p;

} quadElement;


#endif