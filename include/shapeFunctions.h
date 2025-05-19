#ifndef SHAPE_FUNCTIONS_H
#define SHAPE_FUNCTIONS_H

// 2D Shape Functions
typedef struct
{
    double *Ni;
    double *NiPxi;
    double *NiPeta;
    double weight;
} shapeFunctions2D;


int quadShapeFunction (shapeFunctions2D *sf);
int triangularShapeFunction ();
void print2DShapeFunctions (shapeFunctions2D *sf, int gp);

// 3D Shape Functions
int hexShapeFunction ();
int tetraShapeFunction ();

#endif