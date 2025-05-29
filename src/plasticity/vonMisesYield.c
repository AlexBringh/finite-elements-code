#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "vonMisesYield.h"
#include "matrixArithmetic.h"


int trialStress (double *sigmaTrial, double *De, double *epsilon, double *epsilonP, int Dn, int gpCurrent)
{
    /*
        Calculates the trial stress at a Gauss Point adjusted for the accumulated plastic strain tensor.

        Inputs:
        double *sigmaTrial -> Pointer to storing the trial stress vector
        double *De         -> Pointer to D matrix for Hooke's material law (elastic prediction)
        double *epsilon    -> Pointer to strain vector calculated from displacement and B matrix for the current Gauss Point
        double *epsilonP   -> Pointer to plastic strain vector accumulated for all the Gauss Points in the element
        int Dn             -> size of De matrix.
    */

    double *epsilonCorrected = malloc(Dn * sizeof(double));
    double *epsilonPGP = malloc(Dn * sizeof(double));

    // Reset trial stress
    for (int i = 0; i < Dn * Dn; i++) *(sigmaTrial + i) = 0;

    // Extract the plsatic strain tensors for the current Gauss Point
    for (int i = 0; i < Dn; i++)
    {
        *(epsilonCorrected + i) = 0.0;
        *(epsilonPGP + i) = *(epsilonP + gpCurrent * Dn + i);
    }

    // Calculate epsilon - epsilon_p
    matrixSubtract(epsilon, epsilonPGP, epsilonCorrected, Dn, 1);

    // Calculate De * (epsilon - epsilon_p)
    matrixMultiply(De, epsilonCorrected, sigmaTrial, Dn, Dn, 1);

    free(epsilonPGP);
    free(epsilonCorrected);
}

int deviatoricStress2D (double *sDeviatoric, double *sigmaTrial)
{
    /*
        Calculates the deviatoric stress and stores it.

        Inputs:
        double *sDeviatoric -> Pointer to deviatoric stress vector. (Must be 3x1 vector)
        double *sigmaTrial  -> Pointer to trial stress vector. (Must be 3x1 vector)
    */

    *(sDeviatoric + 0) = *(sigmaTrial + 0) - (*(sigmaTrial + 0) + *(sigmaTrial + 1)) / 3;
    *(sDeviatoric + 1) = *(sigmaTrial + 1) - (*(sigmaTrial + 0) + *(sigmaTrial + 1)) / 3;
    *(sDeviatoric + 2) = *(sigmaTrial + 2);
}

int unitDeviatoricStress (double *nDeviatoric, double *sDeviatoric, int sn)
{
    /*
        Calculates the n-unit deviatoric stress and stores it in given pointer.

        Inputs:
        double *nDeviatoric -> Pointer to n-unit deviatoric stress. (Must be (sn)x1 vector)
        double *sDeviatoric -> Pointer to s-deviatoric stress. (Must be (sn)x1 vector)
    */

    double sAbs = 0.0;

    for (int i = 0; i < sn; i++)
    {
        sAbs += pow(*(sDeviatoric + i), 2);
    }
    sAbs = sqrt(sAbs);

    for (int i = 0; i < sn; i++)
    {
        *(nDeviatoric + i) = *(sDeviatoric + i) / sAbs;
    }
}

double vonMisesEquivalentStress2D (double *sDeviatoric)
{
    /*
        Calculates the von Mises equivalent stress based on the s-deviatoric stress

        Input:
        double *sDeviatoric -> Pointer to s-deviatoric stress vector. (Must be a 3x1 vector)

        Output:
        double sigmaEq -> calculated von Mises equivalent stress.
    */
    return sqrt( 3.0 / 2.0 * ( pow(*(sDeviatoric + 0), 2.0) + pow(*(sDeviatoric + 1), 2.0) + 2.0 * pow(*(sDeviatoric + 2), 2.0) ) );
}

double plasticCorrectedYieldStress (double sigmaYieldInitial, double H, double epsilonBarP)
{
    /*
        Calculates the yield stress for the Gauss Point, adjusted for Work Hardening from accumulated equivalent plastic strain.

        Inputs:
        double sigmaYieldInitial -> The initial yield stress for the material at the Gauss Point
        double H                 -> Hardening modulus
        double epsilonBarP       -> equivalent plastic strain accumulated at the Gauss Point

        Output:
        double sigmaYield -> adjusted yield stress at the Gauss Point.
    */

    return sigmaYieldInitial + H * epsilonBarP;
}

double vonMisesYieldFunction (double sigmaEq, double sigmaYield)
{
    /*
        Calculates the value of the von Mises Yield function as sigma_eq - sigma_y.

        Inputs:
        double sigmaEq    -> The von Mises equivalent stress for the current Gauss Point.
        double sigmaYield -> The yield stress, adjusted for work hardening.

        Output:
        double vonMisesYieldFunction -> Returns the yield function value.
    */
    return sigmaEq - sigmaYield;
}