#include <stdio.h>
#include <stdlib.h>

#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"

int SetFieldDensityType(void)
{

  // int FieldDensityType[MAX_NUMBER_OF_FIELD_TYPES];

  FieldDensityType[Density]               = 1;
  FieldDensityType[TotalEnergy]           = 0;
  FieldDensityType[InternalEnergy]        = 0;
  FieldDensityType[Pressure]              = 0;
  FieldDensityType[Velocity1]             = 0;
  FieldDensityType[Velocity2]             = 0;
  FieldDensityType[Velocity3]             = 0;
  FieldDensityType[ElectronDensity]       = 1;
  FieldDensityType[HIDensity]             = 1;
  FieldDensityType[HIIDensity]            = 1;
  FieldDensityType[HeIDensity]            = 1;
  FieldDensityType[HeIIDensity]           = 1;
  FieldDensityType[HeIIIDensity]          = 1;
  FieldDensityType[HMDensity]             = 1;
  FieldDensityType[H2IDensity]            = 1;
  FieldDensityType[H2IIDensity]           = 1;
  FieldDensityType[DIDensity]             = 1;
  FieldDensityType[DIIDensity]            = 1;
  FieldDensityType[HDIDensity]            = 1;
  FieldDensityType[Metallicity]           = 1;
  FieldDensityType[ExtraType0]            = 1;
  FieldDensityType[ExtraType1]            = 1;
  FieldDensityType[GravPotential]         = 0;
  FieldDensityType[Acceleration0]         = 0;
  FieldDensityType[Acceleration1]         = 0;
  FieldDensityType[Acceleration2]         = 0;

  FieldDensityType[CIDensity]             = 1;
  FieldDensityType[CIIDensity]            = 1;
  FieldDensityType[OIDensity]             = 1;
  FieldDensityType[OIIDensity]            = 1;
  FieldDensityType[SiIDensity]            = 1;
  FieldDensityType[SiIIDensity]           = 1;
  FieldDensityType[SiIIIDensity]          = 1;
  FieldDensityType[CHIDensity]            = 1;
  FieldDensityType[CH2IDensity]           = 1;
  FieldDensityType[CH3IIDensity]          = 1;
  FieldDensityType[C2IDensity]            = 1;
  FieldDensityType[COIDensity]            = 1;
  FieldDensityType[HCOIIDensity]          = 1;
  FieldDensityType[OHIDensity]            = 1;
  FieldDensityType[H2OIDensity]           = 1;
  FieldDensityType[O2IDensity]            = 1;

// What is the logic here?
// CRModel can be 0, 1, 2 or 3
// Mach + 2 + CRModel - 1    =    Mach + 1,   Mach + 2, Mach + 3 or Mach + 4  (BAD)
  FieldDensityType[Mach]                  = 0;
  FieldDensityType[PreShockTemperature]   = 0;
  FieldDensityType[PreShockDensity]       = 0;
  FieldDensityType[CRDensity]             = 0;

  FieldDensityType[gParticlePosition]     = 0;
  FieldDensityType[gParticleVelocity]     = 0;
  FieldDensityType[gParticleMass]         = 0;
  FieldDensityType[gParticleAcceleration] = 0;
  FieldDensityType[gParticleNumber]       = 0;
  FieldDensityType[gParticleType]         = 0;
  FieldDensityType[gParticleAttribute]    = 0;
  FieldDensityType[gPotentialField]       = 0;
  FieldDensityType[gAccelerationField]    = 0;
  FieldDensityType[gGravitatingMassField] = 0;
  FieldDensityType[gFlaggingField]        = 0;
  FieldDensityType[gVelocity]             = 0;

  FieldDensityType[RadiationFreq0]        = 0;
  FieldDensityType[RadiationFreq1]        = 0;
  FieldDensityType[RadiationFreq2]        = 0;
  FieldDensityType[RadiationFreq3]        = 0;
  FieldDensityType[RadiationFreq4]        = 0;
  FieldDensityType[RadiationFreq5]        = 0;
  FieldDensityType[RadiationFreq6]        = 0;
  FieldDensityType[RadiationFreq7]        = 0;
  FieldDensityType[RadiationFreq8]        = 0;
  FieldDensityType[RadiationFreq9]        = 0;

  FieldDensityType[EmissivityField0]      = 0;
  FieldDensityType[EmissivityField1]      = 0;

  FieldDensityType[kphHI]                 = 1;
  FieldDensityType[PhotoGamma]            = 0;
  FieldDensityType[kphHeI]                = 1;
  FieldDensityType[kphHeII]               = 1;
  FieldDensityType[kdissH2I]              = 1;

  FieldDensityType[FieldUndefined]        = 0;

  return SUCCESS;

}
