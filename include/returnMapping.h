#ifndef RETURN_MAPPING_H
#define RETURN_MAPPING_H

int elastoPlasticDMatrix (double* D, double* n, double H, int Dn);
double shearModulus (double E, double v);
double plasticMultiplier (double f, double G, double H);

#endif