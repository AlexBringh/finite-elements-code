#include <math.h>

#include "shapeFunctions.h"

int quadShapeFunction (shapeFunctions2D *sf)
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


    // Define constants for 2D quad elements.
    int gp = 4;
    int DOF = 2;
    int dim = 2;
    int nnodes = 4;

    // Allocate memory for the struct
    for (int i = 0; i < gp; i++)
    {
        sf[i].Ni = malloc(nnodes * sizeof(double));
        sf[i].NiPxi = malloc(nnodes * sizeof(double));
        sf[i].NiPeta = malloc(nnodes * sizeof(double));
    }

    double *xi = malloc(gp * sizeof(double));
    double *eta = malloc(gp * sizeof(double));

    *(xi + 0) = -1 / sqrt(3);
    *(xi + 1) =  1 / sqrt(3);
    *(xi + 2) =  1 / sqrt(3);
    *(xi + 3) = -1 / sqrt(3);

    *(eta + 0) = -1 / sqrt(3);
    *(eta + 1) = -1 / sqrt(3);
    *(eta + 2) =  1 / sqrt(3);
    *(eta + 3) =  1 / sqrt(3);


    // Loop through each Gauss Point
    for (int i = 0; i < gp; i++)
    {
        sf[i].Ni[0] = 1 / 4 * (1 - *(xi + i)) * (1 - *(eta + i)); // N1
        sf[i].Ni[1] = 1 / 4 * (1 + *(xi + i)) * (1 - *(eta + i)); // N2
        sf[i].Ni[2] = 1 / 4 * (1 + *(xi + i)) * (1 + *(eta + i)); // N3
        sf[i].Ni[3] = 1 / 4 * (1 - *(xi + i)) * (1 + *(eta + i)); // N4

        sf[i].NiPxi[0] = - 1 / 4 * (1 - *(eta + i)); // N1Pxi
        sf[i].NiPxi[1] =   1 / 4 * (1 - *(eta + i)); // N2Pxi
        sf[i].NiPxi[2] =   1 / 4 * (1 + *(eta + i)); // N3Pxi
        sf[i].NiPxi[3] = - 1 / 4 * (1 + *(eta + i)); // N4Pxi

        sf[i].NiPeta[0] = - 1 / 4 * (1 - *(xi + i)); // N1Peta
        sf[i].NiPeta[1] = - 1 / 4 * (1 + *(xi + i)); // N2Peta
        sf[i].NiPeta[2] =   1 / 4 * (1 + *(xi + i)); // N3Peta
        sf[i].NiPeta[3] =   1 / 4 * (1 - *(xi + i)); // N4Peta

        sf[i].weight = 1.0; // w, integration weight
    }

    free(xi);
    free(eta);
}

int triangularShapeFunction ()
{
    /*
        Not yet implemented for this code. Only 2D quad elements are supported for now.
    */
}

void print2DShapeFunctions (shapeFunctions2D *sf, int gp)
{
    printf("Shape functions for %1d Gauss Point(s). \nNi\tNiPxi\tNiPeta\tWeight\n", gp);
    for (int i = 0; i < gp; i++)
    {
        //double NiValue = *(Ni + i);
        //double NiPxiValue = *(NiPxi + i);
        //double NiPetaValue = *(NiPeta + i);
        //double weightValue = *(weights + i);
        //printf("%.4f\t%.4f\t%.4f\t%.4f\n");
    }
    printf("End of shape functions. \n\n");
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
