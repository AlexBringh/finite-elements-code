#include "inputData.h"
#include "elements.h"

int inputData (quadElement *element, int nElements, int nnodesElement)
{
    int dof = 2;
    int *nodeids;
    double *coords;

    for (int i = 0; i < nElements; i++)
    {
        initQuadElement(&element[i], i, 4, 4, dof, nodeids, coords);
    }
}