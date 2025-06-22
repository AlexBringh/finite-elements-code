#include "elastic.h"

/*

    This module implements functions for calculating elastic stress and strain at a Gauss Point.
    It does not include the logic behind assembling the element- or global stiffness matrices, or internal force vectors.
    
    See module src/elements/elements.c      -> element stiffness matrix and element internal force vector.
    See module src/system/stiffnessMatrix.c -> global stiffness matrix
    See module src/system/forceVector.c     -> global internal force vector

*/

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