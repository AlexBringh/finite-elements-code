#include "elements.h"

/*

    This module implements any functions for initializing and modifying element properties inside the structs defined in the header file.

    See include/elements.h  -> Element struct definition

*/

char* elementType (element e)
{
    /*
        Checks the parsed element struct and determines what element type it is based on the dimensionality and number of nodes.

        Input:
        element e   -> Pointer to element struct

        Output:
        char* type  -> String with the name of the element.
    */

    // TODO: Test if this way of returning char* is OK.

    // 1D elements
    if (e.dim == 1)
    {
        if (e.nnodes == 1) return "1D Point element \n";
        else if(e.nnodes > 1) return "1D Line element \n";
    }
    
    // 2D elements
    else if (e.dim == 2)
    {
        if (e.nnodes == 3) return "2D Triangle element \n";
        else if (e.nnodes == 4) return "2D Quadrilateral element \n";
        else if (e.nnodes == 6) return "2D Quadratic Triangle element \n";
        else if (e.nnodes == 8) return "2D Quadratic Quad element \n";
        else if (e.nnodes == 9) return "2D Full Quadratic element \n";
        else if (e.nnodes == 10) return "2D Cubic Triangle element \n";
    }

    // 3D elements
    else if (e.dim == 3)
    {
        if (e.nnodes == 4) return "3D Linear Tetraherdon element \n";
        else if (e.nnodes == 8) return "3D Linear Hexahedron element \n";
        else if (e.nnodes == 10) return "3D Quadratic Tetrahedron element \n";
        else if (e.nnodes == 20) return "3D Quadratic Hexahedron element \n";
        else if (e.nnodes == 27) return "3D Full Quadratic Hexahedron element \n";
    }


    return NULL; // Undefined
}

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