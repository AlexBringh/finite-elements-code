#include <stdlib.h>

#include "matrixArithmetic.h"
#include "forceVector.h"

void initElementInternalForceVector (double *FintE, int Kem)
{
    /*
        Sets all values of the element internal force vector to 0. Used for initializing and clearing the element internal force vector.
        Cannot be used for Skylie Matrix structures.

        Inputs:
        double *FintE -> Pointer to element internal force vector in array-form.
        int Kem       -> Number of columns (rows) of the element stiffness matrix.
    */
    for (int i = 0; i < Kem; i++)
    {
        *(FintE + i) = 0;
    }
}

void elementInternalForceVector (double *FintE, double *Btrans, double *sigma, double detJ, double weight, int Kem, int Bn)
{
    /*
        Calculates the element internal force vector contribution for a given Gauss Point.
        This function will be ran once for every Gauss Point of an element, with its own unique B transpose, detJ, and weight (though integration weight is not always unique. See 2D quad elements.)
        The function adds the contribution to the element internal force vector for the current element. The vector should be reset at the start of every element in the system, but not for every Gauss Point.

        Input: 
        double *FintE  -> Pointer to element internal force vector. Results are stored here
        double *Btrans -> Pointer to B transpose matrix
        double detJ    -> Determinant value of the Jacobian for the current Gauss Point
        double weight  -> weight value for the current Gauss Point
        int Kem        -> Size of the element stiffness matrix rows. The internal force vector must correspond to this with it's own rows.
        int Bn         -> Size of the columns of the B transpose.
    */

    // Make temporary variable to store the matrix multiplication, Btrans x sigma
    double *Ftemp = malloc(Kem * 1 * sizeof(double));

    // Matrix multiply Btrans with sigma
    matrixMultiply(Btrans, sigma, Ftemp, Kem, Bn, 1);

    // Add the results of Ftemp to FintE, and multiply in the determinant and integration weight.
    for (int i = 0; i < Kem; i++)
    {
        *(FintE + i) += *(Ftemp + i) * detJ * weight;
    }
}


void initGlobalInternalForceVector (double *Fint, int Km)
{
    /*
        Sets all values of the global internal force vector to 0. Used for initializing and clearing the global internal force vector.
        Cannot be used for Skylie Matrix structures.

        Inputs:
        double *Fint -> Pointer to global internal force vector in array-form.
        int Km       -> Number of columns (rows) of the global stiffness matrix.
    */
    for (int i = 0; i < Km; i++)
    {
        *(Fint + i) = 0;
    }
}


void globalInternalForceVector (double *Fint, double *FintE, int *nodeids, int dof, int nnodesElement, int Km, int Kem)
{
    /*
        Assembles the element internal force vector into the global internal force vector. 
        Determines the correct indices based on the node id found in the current element.

        Inputs:
        double *Fint          ->  Pointer to global internal force vector. Results are stored here.
        double *FintE         ->  Pointer to element internal force vector. 
        int *nodeids       ->  Pointer to nodeids for the nodes building the current element.
        int dof            ->  Degrees of freedom of the nodes.
        int nnodesElement  ->  Number of nodes in the element.
        int Km             ->  Size, m_global, of the rows/columns of the global stiffness matrix. The global internal force vector's row indices correspond to this.
        int Kem            ->  Size, m_element, of the rows/columns of the element stiffness matrix. The element internal force vector's row indices correspond to this.
    */

    // Find mapping between element indices and gloval indices
    int *globDOFs = malloc (Kem * sizeof(int));
    for (int i = 0; i < nnodesElement; i++)
    {
        for (int d = 0; d < dof; d++)
        {
            *(globDOFs + i * dof + d) = *(nodeids + i) * dof + d;
        }
    }

    int Fi; // Row indices for the global internal force vector

    // Assemble the element- to global internal force vector.
    for (int i = 0; i < Kem; i++)
    {
        Fi = *(globDOFs + i);

        // Fint[Fi] = FintE[i]
        *(Fint + Fi) += *(FintE + i);
    }

    free(globDOFs);
}

void applyFixedDisplacementResidualVector (int *uFixed, double *r, int Km)
{

}
