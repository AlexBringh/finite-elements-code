#include "elements.h"
#include "matrixArithmetic.h"
#include <stdio.h>

void commitTrialValuesAtGaussPoints (quadElement *element, int nelements)
{
    /*
        Moves the trial values to commited values for all Gauss Points in all elements.

        Inputs:
        quadElement *element  ->  Pointer to all element structs. 
        int nelements         ->  Number of elements.
    */

    for (int e = 0; e < nelements; e++)
    {
        for (int i = 0; i < element[e].gp; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                element[e].sigma[i * 3 + j] = element[e].trialSigma[i * 3 + j];
                element[e].epsilonP[i * 3 + j] = element[e].trialEpsilonP[i * 3 + j];
            }
            element[e].epsilonBarP[i] = element[e].trialEpsilonBarP[i];
        }
    }
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