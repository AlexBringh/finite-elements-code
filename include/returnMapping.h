#ifndef RETURN_MAPPING_H
#define RETURN_MAPPING_H

int elastoPlasticDMatrix (double* D, double* n, double H, int Dn);
double shearModulus (double E, double v);
double plasticMultiplier (double f, double G, double H);
void plasticStressCorrection (double *sigma, double *sigmaTrial, double *n, double G, double deltaGamma, int index);
void trialPlasticStrain (double* trialEpsilonP, double *epsilonP, double *n, double deltaGamma, int index, int gpCurrent);
double trialEquivalentPlasticStrain (double deltaGamma);

#endif