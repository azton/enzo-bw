#ifndef __typedefs_h_
#define __typedefs_h_
/***********************************************************************
/
/  MISCELANEOUS TYPEDEFS AND ENUMERATIONS
/
/  written by: Greg Bryan
/  date:       November, 1994
/  modified1:  S.W. Skillman
/  date:       January, 2009 - Added CosmicRayData.h
/  modified2:  Robert Harkness
/  date:       February, 2010 - Added Radhydro fields
/
/  PURPOSE:
/
************************************************************************/

#include "CloudyCoolingData.h"
#include "CoolData.h"
#include "RateData.h"
#include "RadiationFieldData.h"
#include "TestProblemData.h"
#include "CosmicRayData.h"

/* These are the different types of baryon fields. */

#ifdef SMALL_INTS
typedef int field_type;
typedef int boundary_type;
typedef int gravity_boundary_type;
typedef int interpolation_type;
typedef int hydro_method;
typedef int enum_type;
#endif

#ifdef LARGE_INTS
typedef long_int field_type;
typedef long_int boundary_type;
typedef long_int gravity_boundary_type;
typedef long_int interpolation_type;
typedef long_int hydro_method;
typedef long_int enum_type;
#endif

const field_type 
  Density         = 0,
  TotalEnergy     = 1,
  InternalEnergy  = 2,
  Pressure        = 3,
  Velocity1       = 4,
  Velocity2       = 5,
  Velocity3       = 6,
  ElectronDensity = 7,
  HIDensity       = 8,
  HIIDensity      = 9,
  HeIDensity      = 10,
  HeIIDensity     = 11,
  HeIIIDensity    = 12,
  HMDensity       = 13,
  H2IDensity      = 14,
  H2IIDensity     = 15,
  DIDensity       = 16,
  DIIDensity      = 17,
  HDIDensity      = 18,
  Metallicity     = 19,
  ExtraType0      = 20,
  ExtraType1      = 21,
  GravPotential   = 22,
  Acceleration0   = 23,
  Acceleration1   = 24,
  Acceleration2   = 25,

/* these pseudo-fields are used to access grid data 
   the "g" prefix is to avoid namespace conflict */

  gParticlePosition     = 26,
  gParticleVelocity     = 27,
  gParticleMass         = 28,
  gParticleAcceleration = 29,
  gParticleNumber       = 30,
  gParticleType         = 31,
  gParticleAttribute    = 32,
  gPotentialField       = 33,
  gAccelerationField    = 34,
  gGravitatingMassField = 35,
  gFlaggingField        = 36,
  gVelocity             = 37,

/* these are required for Simon Glover's chemistry (which also needs some of the
   other fields, which are used for MultiSpecies) */
  CIDensity       = 38, 
  CIIDensity      = 39, 
  OIDensity       = 40, 
  OIIDensity      = 41,
  SiIDensity      = 42,
  SiIIDensity     = 43,
  SiIIIDensity    = 44,
  CHIDensity      = 45,
  CH2IDensity     = 46,
  CH3IIDensity    = 47,
  C2IDensity      = 48,
  COIDensity      = 49,
  HCOIIDensity    = 50,
  OHIDensity      = 51,
  H2OIDensity     = 52,
  O2IDensity      = 53,

/* these are required for Sam Skillmans Shock/Cosmic ray models. */
  Mach            = 54,
  PreShockTemperature = 55,
  PreShockDensity = 56,  
  CRDensity       = 57,

/* RT module stuff */
  RadiationFreq0  = 58,
  RadiationFreq1  = 59,
  RadiationFreq2  = 60,
  RadiationFreq3  = 61,
  RadiationFreq4  = 62,
  RadiationFreq5  = 63, 
  RadiationFreq6  = 64,
  RadiationFreq7  = 65, 
  RadiationFreq8  = 66, 
  RadiationFreq9  = 67, 
  
/* Star Maker emissivity fields */
  EmissivityField0 = 68, 
  EmissivityField1 = 74, 

/* these fields are used to provide radiation-induced 
   photo-ionization and photo-heating rates to the 
   chemistry/cooling routines */

  kphHI           = 69,
  PhotoGamma      = 70,
  kphHeI          = 71,
  kphHeII         = 72,
  kdissH2I        = 73,

  FieldUndefined  = 75;  // this needs to be the last one

   
/*
enum field_type {Density, TotalEnergy, InternalEnergy, Pressure,
		 Velocity1, Velocity2, Velocity3, 
		 ElectronDensity, HIDensity, HIIDensity,  HeIDensity, 
		 HeIIDensity, HeIIIDensity, HMDensity, H2IDensity, 
		 H2IIDensity, DIDensity, DIIDensity, HDIDensity,
                 Metallicity, ExtraType0, ExtraType1, GravPotential,
		 Acceleration0, Acceleration1,Acceleration2,
                 FieldUndefined};
*/

//define FieldTypeIsDensity(A) (((A) >= TotalEnergy && (A) <= Velocity3) ? FALSE : TRUE)

/*Define a new function to set either (Mach) or (Mach and CRDensity) to be FALSE.
  This Depends on the way CRModel is defined.
  So check this if you change what each value of CRModel does! */

//define FieldTypeIsDensity(A) ( (  ( (A) >= TotalEnergy && (A) <= Velocity3 ) || ( (A) >= Mach && ( (A) <= (Mach + 2 + CRModel - 1) && ( (A) <= Mach + 3 ) ) ) ) ? FALSE : TRUE)

#define FieldTypeIsDensity(A) (FieldDensityType[(A)])
#define FieldTypeIsColor(A) (FieldColorType[(A)])

/* These are the different types of fluid boundary conditions. */

const boundary_type
  reflecting        = 0,
  outflow           = 1,
  inflow            = 2,
  periodic          = 3,
  BoundaryUndefined = 4;

// enum boundary_type {reflecting, outflow, inflow, periodic, BoundaryUndefined};

/* These are the different types of gravity boundary conditions. */

const gravity_boundary_type
  TopGridPeriodic  = 0,
  TopGridIsolated  = 1,
  SubGridIsolated  = 2,
  GravityUndefined = 3;

// enum gravity_boundary_type {TopGridPeriodic, TopGridIsolated, 
// 				    SubGridIsolated, GravityUndefined};

/* Interpolation types. */

const interpolation_type
  ThirdOrderA            = 0,
  SecondOrderA           = 1,
  SecondOrderB           = 2,
  SecondOrderC           = 3,
  FirstOrderA            = 4,
  InterpolationUndefined = 5;


// enum interpolation_type {ThirdOrderA, SecondOrderA, SecondOrderB, SecondOrderC,
// 			 FirstOrderA, InterpolationUndefined};

/* Hydrodynamics methods. */

const hydro_method
  PPM_DirectEuler      = 0,
  PPM_LagrangeRemap    = 1,
  Zeus_Hydro           = 2,
  HydroMethodUndefined = 3;

// enum hydro_method {PPM_DirectEuler, PPM_LagrangeRemap, Zeus_Hydro};

const enum_type
  PLM                  = 0,
  PPM                  = 1;

const enum_type
  FluxReconstruction   = 0,
  HLL                  = 1,
  Marquina             = 2,
  LLF                  = 3,
  HLLC                 = 4,
  TwoShock             = 5;



/* Define a float/int union. */

union float_int {
  long_int ival;
  float fval;
  FLOAT FVAL;
};

#endif
