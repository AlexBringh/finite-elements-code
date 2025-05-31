#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "postProcessing.h"

void postProcessingElastoPlastic (double *nodalStress, double *nodalStrain, double *nodalPlasticStrain, int nnodes, double *u, quadElement *elements, int nelements, shapeFunctions2D *sf, int m)
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

    printf("Post-Processing #1\n");
    nodalStress = calloc(nnodes, sizeof(double));
    nodalStrain = calloc(nnodes ,sizeof(double));
    nodalPlasticStrain = calloc(nnodes, sizeof(double));
    int stressDim = 3;
    double sigmaEq;
    double epsilonEq;
    int gp = elements->gp;

    // These are temp values
    double *stress = calloc(stressDim * nnodes, sizeof(double));
    double *strain = calloc(stressDim * nnodes, sizeof(double));
    double *weights = calloc(nnodes, sizeof(double));
    double *tempStress = calloc(stressDim, sizeof(double));
    double *tempStrain = calloc(stressDim, sizeof(double));
    double *ue = calloc(elements->nnodes * elements->dof, sizeof(double));
    double *J = calloc(2 * 2, sizeof(double));
    double *Jinv = calloc(2 * 2, sizeof(double));
    double detJ = 0;
    double *B = calloc(stressDim * 2 * elements->nnodes, sizeof(double));
    double *Btrans = calloc(stressDim * 2 * elements->nnodes, sizeof(double));
    double *Ni = calloc(4 * gp, sizeof(double));
    double N;
    int nodeid;
    printf("Post-Processing #2\n");

    // Store the shape functions for convenience and code clarity.
    for (int p = 0; p < gp; p++) 
    {
        for (int i = 0; i < 4; i++)
        {
            *(Ni + p * gp + i) = sf[p].Ni[i];
        }
    }
    printf("Post-Processing #3\n");

    // Interpolate the contribution of each Gauss Point in each element to the nodal stresses and strains.
    for (int e = 0; e < nelements; e++) // Loop over elements
    {
        for (int i = 0; i < elements[e].nnodes; i++)
        {
            if (elements[e].nodeids == NULL) {
                printf("Error: elements[%d].nodeids is NULL\n", e);
                exit(1);
            }
            int node_id = elements[e].nodeids[i];
            if (node_id < 0 || node_id >= nnodes) {
                printf("Error: node_id %d out of bounds for element %d (i=%d)\n", node_id, e, i);
                exit(1);
            }

            for (int d = 0; d < elements[e].dof; d++) // Loop over each DOF
            {
                *(ue + i * elements[e].dof + d) = *(u + elements[e].nodeids[i] * elements[e].dof + d); // Store the displacement from the node with this id.
            }
        }

        for (int i = 0; i < elements[e].nnodes; i++) // Loop over nodes in the element
        {
            nodeid = elements[e].nodeids[i]; // Current node id
            for (int p = 0; p < gp; p++) // Loop over Gauss Points in the element
            {
                // Calculate the strain from the current displacement, since we do not store strain normally.
                quadJacobian(J, Jinv, 2, &detJ, sf[p].Ni, sf[p].NiPxi, sf[p].NiPeta, elements[e].coords, elements[e].dof);
                quadBMatrix(B, Btrans, Jinv, sf[p].NiPxi, sf[p].NiPeta);
                displacementStrain(tempStrain, B, ue, elements[e].nnodes, elements[e].dof);

                // Store the equivalent plastic strain interpolated with the shape functions
                *(nodalPlasticStrain + nodeid) += *(Ni + p * gp + i) * elements[e].epsilonBarP[p];
                *(weights + nodeid) += 1.0; // Add one for every time a Gauss Points contributes to the interpolation. Used for averaging later.

                for (int j = 0; j < stressDim; j++) // Loop over every stress dim. For 2D this is xx (normal), yy (normal) and xy (shear)
                {
                    double NiCurrent = *(Ni + p * gp + i);
                    printf("%.4f \n", elements[e].sigma[p * stressDim + j]);
                    *(stress + nodeid * stressDim + j) += NiCurrent * elements[e].sigma[p * stressDim + j]; // Shape function of every Gauss Point, corresponding to the i-th element, multiplied by the stress at the Gauss Point, for the current dim.
                    *(strain + nodeid * stressDim + j) += NiCurrent * *(tempStrain + j); // Interpolated strain
                    
                }
            }
        }
    }
    printf("Post-Processing #4\n");

    // Average out the contibutions for each node.
    for (int i = 0; i < nnodes; i++)
    {
        if (*(weights + i) > 0)
        {
            for (int j = 0; j < stressDim; j++)
            {
                *(stress + i * stressDim + j) /= (double) *(weights + i);
                *(strain + i * stressDim + j) /= (double) *(weights + i);
            }
            *(nodalPlasticStrain + i) /= (double) *(weights + i);
        }
    }
    printf("Post-Processing #5\n");


    // Convert stress and strain at the nodes to von Mises equivalent stress and strain.
    for (int i = 0; i < nnodes; i++)
    {
        // Store each stress/strain dim in tempStress and tempStrain.
        for (int j = 0; j < stressDim; j++) 
        {
            *(tempStress + j) = *(stress + i * stressDim + j);
            *(tempStrain + j) = *(strain + i * stressDim + j);
            printf("Post-Processing #6: %.4f\n", *(strain + i * stressDim + j));
        }
        
        // Calculate and store von Mises Equivalent stress and strain.
        *(nodalStress + i) = vonMisesEquivalentStress2D(tempStress);
        *(nodalStrain + i) = vonMisesEquivalentStrain2D(tempStrain);
        printf("Post-Processing #7: %.2f\n", *(nodalStress + i));
        printf("Post-Processing #7.5: %.2f\n", *(nodalStrain + i));
    }   
    printf("Post-Processing #8\n");

    free(stress);
    free(strain);
    free(weights);
    free(tempStress);
    free(tempStrain);
    free(ue);
    free(J);
    free(Jinv);
    free(B);
    free(Btrans);
    free(Ni);
}
