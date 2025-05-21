#include "elements.h"
#include "matrixArithmetic.h"

int initQuadElement (quadElement *element, int id, int gp, int nnodesElement, int dof, int *nodeids, double *coords)
{
    /*
        Initializes the struct element, allocates memory for all the internal arrays and stores the internal values. Sets accumulated- and trial strains and stresses to 0.

        Inputs:
        quadElement *element -> Pointer to current element struct. Result is stored here.
        int id               -> Id of the element.
        int gp               -> Number of Gauss Points in the element.
        int nnodesElement    -> Number of nodes in the element.
        int dof              -> Degrees of freedom of the nodes.
        int *nodeids         -> Pointer to node ids of the element nodes.
        int *coords          -> Pointer to global coordinates of the nodes.
    */

    element->id = id;
    element->gp = gp;
    element->nnodes = nnodesElement;
    element->dof = dof;

    
    

    element->nodeids = malloc(nnodesElement * sizeof(int));
    element->coords = malloc(nnodesElement * dof * sizeof(double));
    element->sigma = malloc(3 * gp * sizeof(double));
    element->epsilonP = malloc(3 * gp * sizeof(double));
    element->epsilonBarP = malloc(gp * sizeof(double));
    element->trialSigma = malloc(3 * gp * sizeof(double));
    element->trialEpsilonP = malloc(3 * gp * sizeof(double));
    element->trialEpsilonBarP = malloc(gp * sizeof(double));

    // Store element ids.
    for (int i = 0; i < nnodesElement; i++)
    {
        element->nodeids[i] = *(nodeids + i);

        // Store coords of the nodes
        for (int d = 0; d < dof; d++)
        {
            element->coords[i * dof + d] = *(coords + i * dof + d);
        }
    }

    // Init stress and strain. Each vector has 3 values for each Gauss Point.
    for (int i = 0; i < 3 * gp; i++)
    {
        element->sigma[i] = 0;
        element->trialSigma[i] = 0;
        element->epsilonP[i] = 0;
        element->trialEpsilonP[i] = 0;
    }

    // Init equivalent strain.
    for (int i = 0; i < gp; i++)
    {
        element->epsilonBarP[i] = 0;
        element->trialEpsilonBarP[i] = 0;
    }

    return 0;
}


int displacementStrain (double *epsilon, double *B, double *ue, int nnodesElement, int dof)
{
    /*
        Calculates the strain at the Gauss Point determined by the displacement of the nodes in the element, and the B matrix calculated for the current Gauss Point.

        Inputs:
        double *epsilon   -> Pointer to strain at the current Gauss Point. Results are stored here.
        double *B         -> Pointer to B matrix at the current Gauss Point.
        double *ue        -> Pointer to element displacement vector.
        int nnodesElement -> Number of nodes in the element.
        int dof           -> Degrees of freedom of the nodes.
    */

    // Init / reset the strain tensor, epsilon
    for (int i = 0; i < nnodesElement * dof; i++)
    {
        *(epsilon + i) = 0;
    }

    // Perform the matrix multiplication
    matrixMultiply(B, ue, epsilon, 3, (nnodesElement * dof), 1);

    return 0;
}