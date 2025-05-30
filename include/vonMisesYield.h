#ifndef VON_MISES_YIELD_H
#define VON_MISES_YIELD_H

int trialStress (double *sigmaTrial, double *De, double *epsilon, double *epsilonP, int Dn, int gpCurrent);
int deviatoricStress2D (double *sDeviatoric, double *sigmaTrial);
int unitDeviatoricStress (double *nDeviatoric, double *sDeviatoric, int sn);
double vonMisesEquivalentStress2D (double *sDeviatoric);
double vonMisesEquivalentStrain2D (double *epsilon);
double plasticCorrectedYieldStress (double sigmaYieldInitial, double H, double epsilonBarP);
double vonMisesYieldFunction (double sigmaEq, double sigmaYield);

#endif