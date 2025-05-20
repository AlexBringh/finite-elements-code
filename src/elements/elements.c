#include "elements.h"
#include "matrixArithmetic.h"

int initQuadElement ()
{
    /*
        TODO: Documentation
    */

    return 0;
}


int displacementStrain (double *epsilon, double *B, double *ue, int nnodesElement, int dof)
{
    /*
        TODO: Documentation
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