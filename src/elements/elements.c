#include "elements.h"

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


int elementInit (element e, int gp, int dof, int dim, int nnodes, int mat, int flag)
{
    /*
    IMPLEMENT: Init function for element struct. Uses flag to determine if the simulation will accumulate plastic effects.
    
    */

    e.gp = gp;
    e.dof = dof;
    e.dim = dim;
    e.nnodes = nnodes;
    e.mat = mat;

    if (flag == 1) // Plastic variables
    {

    }

    return 0;
}

int elementSetElasticTrialValues(element e, int gp, double sigma)
{
    /*
        IMPLEMENT: Set trial values for stress
    */

    return 0;
}

int elementSetPlasticTrialValues (element e, int gp, double epsilonP, double epsilonBarP)
{
    /*
        IMPLEMENT: Set trial values for plastic strains
    */

    return 0;
}

void elementCommitElasticTrialValues (element e)
{
    /*
        IMPLEMENT: Commit elastic trial values to accumulated values
    */


}

void elementCommitPlasticTrialValues (element e)
{
    /*
        IMPLEMENT: Commit plastic trial values to accumulated values
    */



}