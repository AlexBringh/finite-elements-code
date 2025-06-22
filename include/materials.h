#ifndef MATERIALS_H
#define MATERIALS_H

typedef struct 
{
    // General material type
    int id;     // Material id
    double E;   // Young' Modulus
    double v;   // Poisson's ratio
    double H;   // Hardening Modulus for Strain/Work Hardening
    double y;   // Material yield stress

} material;

typedef struct 
{
    // Solid material type
    int id;     // Material id
    double E;   // Young' Modulus
    double v;   // Poisson's ratio
    double H;   // Hardening Modulus for Strain/Work Hardening
    double y;   // Material yield stress

} solidMaterial;

typedef struct
{
    // Fluid material type
    int id;
    // IMPLEMENT: fluid material properties

} fluidMaterial;

#endif