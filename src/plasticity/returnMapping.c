#include "returnMapping.h"

int elastoPlasticDMatrix (double* D, double* n, double H, int Dn)
{
    /*
        Computesthe elasto-plastic material stiffness matrix, Dep = De - (De * n)(De * n)^T / (n^T * De * n + 2H/3)
        
        Inputs:
        double *D   -> D, stiffness matrix [Dn x Dn]
        double *n   -> unit deviatoric flow direction vector [Dn]
        double H    -> hardening modulus
        int Dn      -> size of the system (3 for 2D, 6 for 3D)   
    */


    double *De = malloc(Dn * Dn * sizeof(double));
    double *v = malloc(Dn * sizeof(double));
    double *numerator = malloc(Dn * Dn * sizeof(double));
    double Hprime = (2.0 / 3.0) * H;
    double denominator = Hprime;

    // Take the elastic material stiffness matrix value and store it in a temporary memory slot.
    for (int i = 0; i < Dn * Dn; i++) 
    {
        *(De + i) = *(D + i);
    }
    // Calculation begins here.
    // Step 1: Compute v = De * n
    for (int i = 0; i < Dn; i++) {
        v[i] = 0.0;
        for (int j = 0; j < Dn; j++) {
            v[i] += De[i * Dn + j] * n[j];
        }
    }

    // Step 2: Compute numerator = v * v^T
    for (int i = 0; i < Dn; i++) {
        for (int j = 0; j < Dn; j++) {
            numerator[i * Dn + j] = v[i] * v[j];
        }
    }

    // Step 3: Compute denominator = n^T * De * n + H'
    for (int i = 0; i < Dn; i++) {
        denominator += n[i] * v[i];  // since v = De * n
    }

    // Step 4: Compute Dep = De - numerator / denominator
    for (int i = 0; i < Dn; i++) {
        for (int j = 0; j < Dn; j++) {
            D[i * Dn + j] = De[i * Dn + j] - numerator[i * Dn + j] / denominator;
        }
    }

    free(De);
    free(v);
    free(numerator);
}


double shearModulus (double E, double v)
{
    /*
        Calculates the shear modulus, G.

        Inputs:
        double E -> Young's Modulus
        double v -> Poisson's Ratio

        Output:
        double G -> Shear Modulus
    */

    return E / (2 * (1+v));
}


double plasticMultiplier (double f, double G, double H)
{
    /*
        Calculates the plastic multiplier, delta-gamma.

        Inputs:
        double f -> von Mises yield function value
        double G -> Shear Modulus
        double H -> Hardening Modulus

        Output:
        double deltaGamma -> Calculated Plastic Multiplier
    */

    return f / (3 * G + H);
}

int plasticStressCorrection ()  
{

}


int trialPlasticStrain ()
{

}


double trialEquivalentPlasticStrain ()
{

}