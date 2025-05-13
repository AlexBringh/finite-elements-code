#include <math.h>

#include "shapeFunctions.h"

int quadShapeFunction (double *Ni, double *NiPxi, double *NiPeta, double *weights)
{
    /*
        Input value *Ni should be a pointer with allocated enough memory for 4 double values.
        Input value *NiPxi should be a pointer with allocated enough memory for 4 double values.
        Input value *NiPeta should be a pointer with allocated enough memory for 4 double values.

        This doesn't look so pretty and fancy since the calculations are so similar, but have slightly varying values.
        However some efficiency and beauty is sacrificed for the sake of clarity of where the values are comming from.
        The calculations of the local shape functions is only done once, anyway.
    */
    // TODO: Documentation

    // Define Gaussian Weights for quad elements in 2D space.
    int gp = 4;
    int DOF = 2;
    
    // Weights
    for (int i = 0; i < gp; i++)
    {
        *(weights + i) = 1.0; // All four Gauss points will have the value of 1.0.
    }

    // Gauss Points
    double xi1  =  -1 / sqrt(3); // xi of the first Gauss point
    double eta1 =  -1 / sqrt(3); // eta of the first Gauss point
    double xi2  =   1 / sqrt(3); // xi of the second Gauss point
    double eta2 =  -1 / sqrt(3); // eta of the second Gauss point
    double xi3  =   1 / sqrt(3); // xi of the third Gauss point
    double eta3 =   1 / sqrt(3); // eta of the third Gauss point
    double xi4  =   1 / sqrt(3); // xi of the fourth Gauss point
    double eta4 =  -1 / sqrt(3); // eta of the fourth Gauss point

    // Calculate the local shape functions of a 2D quad mesh element.
    *(Ni + 0) = 1 / 4 * (1 - xi1) * (1 - eta1); // N1
    *(Ni + 1) = 1 / 4 * (1 + xi2) * (1 - eta2); // N2
    *(Ni + 2) = 1 / 4 * (1 + xi3) * (1 + eta3); // N3
    *(Ni + 3) = 1 / 4 * (1 - xi4) * (1 - eta4); // N4

    // Calculate the partial derivatives of the local shape functions with respects to xi.
    *(NiPxi + 0) = -1 / 4 * (1 - eta1); // N1Pxi
    *(NiPxi + 1) =  1 / 4 * (1 - eta2); // N2Pxi
    *(NiPxi + 2) =  1 / 4 * (1 + eta3); // N3Pxi
    *(NiPxi + 3) = -1 / 4 * (1 + eta4); // N4Pxi

    // Calculate the partial derivatives of the local shape functions with respects to eta.
    *(NiPeta + 0) = -1 / 4 * (1 - xi1); // N1Peta
    *(NiPeta + 1) = -1 / 4 * (1 + xi2); // N1Peta
    *(NiPeta + 2) =  1 / 4 * (1 + xi3); // N1Peta
    *(NiPeta + 3) =  1 / 4 * (1 - xi4); // N1Peta


    free(weights);
}

int triangularShapeFunction ()
{
    /*
        Not yet implemented for this code. Only 2D quad elements are supported for now.
    */
}

int hexShapeFunction ()
{
    /*
        Not yet implemented for this code. Only 2D quad elements are supported for now.
    */
}

int tetraShapeFunction ()
{
    /*
        Not yet implemented for this code. Only 2D quad elements are supported for now.
    */
}