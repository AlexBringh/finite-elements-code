#include "jacobian.h"
#include "shapeFunctions.h"


int quadJacobian (double *J, double *Jinv, int m, double *detJ, double *Ni, double *NiPxi, double *NiPeta, double *xi, double *yi)
{
    /*
        TODO: Documentation.
    */
    if (m != 2) return -1; // If the Jacobian is not a 2x2 matrix, then this function should not run.

    // Empty / initialize the Jacobian matrix to avoid problems.
    for (int i = 0; i < m * m; i++)
    {
        *(J + i) = 0;
    }

    // Empty / initialize the Jacobian matrix to avoid problems.
    for (int i = 0; i < m * m; i++)
    {
        *(Jinv + i) = 0;
    }
    
    // Calculate values and assemble the Jacobian matrix.
    for (int i = 0; i < 4; i++)
    {
        *(J + 0) += *(NiPxi + i) * *(xi + i); // Partial derivative of x with respects to Xi (greek letter).
        *(J + 1) += *(NiPxi + i) * *(yi + i); // Partial derivative of y with respects to Xi (greek letter).
        *(J + 2) += *(NiPeta + i) * *(xi + i); // Partial derivative of x with respects to Eta.
        *(J + 3) += *(NiPeta + i) * *(yi + i); // Partial derivative of y with respects to Eta.
    }

    *detJ = *(J + 0) * *(J + 3) - *(J + 2) * *(J + 1); // Determinant of a 2x2 matrix is J_0 * J_3 - J_2 * J_1

    // Assemble the inverse Jacobian matrix.
    *(Jinv + 0) =  1 / *(detJ) * *(J + 3); 
    *(Jinv + 1) = -1 / *(detJ) * *(J + 1);
    *(Jinv + 2) = -1 / *(detJ) * *(J + 2);
    *(Jinv + 3) =  1 / *(detJ) * *(J + 0); 
    
    return 0;
}