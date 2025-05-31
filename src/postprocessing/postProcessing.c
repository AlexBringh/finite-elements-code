#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "postProcessing.h"

void postProcessingElastoPlastic (double *nodalStress, double *nodalStrain, double *nodalPlasticStrain, int nnodes, double *u, quadElement *elements, int nelements, shapeFunctions2D *sf, int m, double E)
{
    /*
        Post-processing for the results. Interpolates the stresses and strains from the Gauss Points to the nodes.

        Inputs:
        double *nodalStress         -> Pointer to nodal von Mises equivalent stress vector. Stores stress for every node (1 stress value per node)
        double *nodalStrain         -> Pointer to nodal von Mises equivalent strain vector (effective strain). Stores strain for every node (1 strain value per node)
        double *nodalPlasticStrain  -> Pointer to nodal equivalent plastic strain vector. Stores equivalent plastic strain for every node (1 value per node)  
        int nnodes                  -> Number of nodes
        double *u                   -> Pointer to global displacement vector for nodal displacement
        quadElement *elements       -> Pointer to elements structs
        int nelements               -> Number of elements
        shapeFunctions2D            -> Pointer to 2D shape function structs
        int m                       -> Size m, number of global degrees of freedom
    */
    int stressDim = 3;
    double sigmaEq;
    double epsilonEq;
    int gp = elements->gp;

    // These are temp values
    double *stress = calloc(stressDim * nnodes, sizeof(double));
    double *weights = calloc(nnodes, sizeof(double));
    double *tempStress = calloc(stressDim, sizeof(double));
    double *J = calloc(2 * 2, sizeof(double));
    double *Jinv = calloc(2 * 2, sizeof(double));
    double detJ = 0;
    double *Ni = calloc(4 * gp, sizeof(double));
    double N;
    int nodeid;

    // Store the shape functions for convenience and code clarity.
    for (int p = 0; p < gp; p++) 
    {
        for (int i = 0; i < 4; i++)
        {
            *(Ni + p * gp + i) = sf[p].Ni[i];
        }
    }

    // Interpolate the contribution of each Gauss Point in each element to the nodal stresses and strains.
    for (int e = 0; e < nelements; e++) // Loop over elements
    {
        for (int i = 0; i < elements[e].nnodes; i++) // Loop over nodes in the element
        {
            nodeid = elements[e].nodeids[i]; // Current node id
            for (int p = 0; p < gp; p++) // Loop over Gauss Points in the element
            {
                // Calculate the strain from the current displacement, since we do not store strain normally.
                quadJacobian(J, Jinv, 2, &detJ, sf[p].Ni, sf[p].NiPxi, sf[p].NiPeta, elements[e].coords, elements[e].dof);

                // Store the equivalent plastic strain interpolated with the shape functions
                *(nodalPlasticStrain + nodeid) += *(Ni + p * gp + i) * elements[e].epsilonBarP[p];
                double NiCurrent = *(Ni + p * gp + i);
                *(weights + nodeid) += sf[p].weight * NiCurrent * detJ; // Add one for every time a Gauss Points contributes to the interpolation. Used for averaging later.

                for (int j = 0; j < stressDim; j++) // Loop over every stress dim. For 2D this is xx (normal), yy (normal) and xy (shear)
                {
                    
                    *(stress + nodeid * stressDim + j) += NiCurrent * elements[e].sigma[p * stressDim + j] * detJ * sf[p].weight; // Shape function of every Gauss Point, corresponding to the i-th element, multiplied by the stress at the Gauss Point, for the current dim.
                }
            }
        }
    }

    for (int i = 0; i < nnodes; i++)
    {
        if (*(weights + i) > 0)
        {
            for (int j = 0; j < stressDim; j++)
            {
                *(stress + i * stressDim + j) /= (double) *(weights + i);
            }
            *(nodalPlasticStrain + i) /= (double) *(weights + i);
        }
    }

    // Convert stress and strain at the nodes to von Mises equivalent stress and strain.
    for (int i = 0; i < nnodes; i++)
    {
        // Store each stress/strain dim in tempStress and tempStrain.
        for (int j = 0; j < stressDim; j++) 
        {
            *(tempStress + j) = *(stress + i * stressDim + j);
        }
        
        // Calculate and store von Mises Equivalent stress and strain.
        *(nodalStress + i) = vonMisesEquivalentStress2D(tempStress);
        *(nodalStrain + i) = *(nodalStress + i) / E;
    }  

    free(stress);
    free(weights);
    free(tempStress);
    free(J);
    free(Jinv);
    free(Ni);
}
