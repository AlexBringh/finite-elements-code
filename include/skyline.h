#ifndef SKYLINE_H
#define SKYLINE_H

#include <stdio.h>
#include <stdlib.h>

skylineMatrix* initSkylineMatrix (int n);
int addSkylineElement (skylineMatrix* matrix, int m, int n, double val);
int getSkylineElement (skylineMatrix* matrix, int m, int n);
int updateSkylineElement (skylineMatrix* matrix, int m, int n, double val);

#endif