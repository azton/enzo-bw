#include <stdio.h>
#include <stdlib.h>

#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"

int SetFieldColorType(void)
{

  // int FieldColorType[MAX_NUMBER_OF_FIELD_TYPES];

  FieldColorType[Density]               = 0;
  FieldColorType[TotalEnergy]           = 0;
  FieldColorType[InternalEnergy]        = 0;
  FieldColorType[Pressure]              = 0;
  FieldColorType[Velocity1]             = 0;
  FieldColorType[Velocity2]             = 0;
  FieldColorType[Velocity3]             = 0;
  FieldColorType[ElectronDensity]       = 1;
  FieldColorType[HIDensity]             = 1;
  FieldColorType[HIIDensity]            = 1;
  FieldColorType[HeIDensity]            = 1;
  FieldColorType[HeIIDensity]           = 1;
  FieldColorType[HeIIIDensity]          = 1;
  FieldColorType[HMDensity]             = 1;
  FieldColorType[H2IDensity]            = 1;
  FieldColorType[H2IIDensity]           = 1;
  FieldColorType[DIDensity]             = 1;
  FieldColorType[DIIDensity]            = 1;
  FieldColorType[HDIDensity]            = 1;
  FieldColorType[Metallicity]           = 1;
  FieldColorType[ExtraType0]            = 1;
  FieldColorType[ExtraType1]            = 1;
  FieldColorType[GravPotential]         = 0;
  FieldColorType[Acceleration0]         = 0;
  FieldColorType[Acceleration1]         = 0;
  FieldColorType[Acceleration2]         = 0;

  FieldColorType[CIDensity]             = 1;
  FieldColorType[CIIDensity]            = 1;
  FieldColorType[OIDensity]             = 1;
  FieldColorType[OIIDensity]            = 1;
  FieldColorType[SiIDensity]            = 1;
  FieldColorType[SiIIDensity]           = 1;
  FieldColorType[SiIIIDensity]          = 1;
  FieldColorType[CHIDensity]            = 1;
  FieldColorType[CH2IDensity]           = 1;
  FieldColorType[CH3IIDensity]          = 1;
  FieldColorType[C2IDensity]            = 1;
  FieldColorType[COIDensity]            = 1;
  FieldColorType[HCOIIDensity]          = 1;
  FieldColorType[OHIDensity]            = 1;
  FieldColorType[H2OIDensity]           = 1;
  FieldColorType[O2IDensity]            = 1;

// What is the logic here?
// CRModel can be 0, 1, 2 or 3
// Mach + 2 + CRModel - 1    =    Mach + 1,   Mach + 2, Mach + 3 or Mach + 4  (BAD)
  FieldColorType[Mach]                  = 0;
  FieldColorType[PreShockTemperature]   = 0;
  FieldColorType[PreShockDensity]       = 0;
  FieldColorType[CRDensity]             = 0;

  FieldColorType[gParticlePosition]     = 0;
  FieldColorType[gParticleVelocity]     = 0;
  FieldColorType[gParticleMass]         = 0;
  FieldColorType[gParticleAcceleration] = 0;
  FieldColorType[gParticleNumber]       = 0;
  FieldColorType[gParticleType]         = 0;
  FieldColorType[gParticleAttribute]    = 0;
  FieldColorType[gPotentialField]       = 0;
  FieldColorType[gAccelerationField]    = 0;
  FieldColorType[gGravitatingMassField] = 0;
  FieldColorType[gFlaggingField]        = 0;
  FieldColorType[gVelocity]             = 0;

  FieldColorType[RadiationFreq0]        = 0;
  FieldColorType[RadiationFreq1]        = 0;
  FieldColorType[RadiationFreq2]        = 0;
  FieldColorType[RadiationFreq3]        = 0;
  FieldColorType[RadiationFreq4]        = 0;
  FieldColorType[RadiationFreq5]        = 0;
  FieldColorType[RadiationFreq6]        = 0;
  FieldColorType[RadiationFreq7]        = 0;
  FieldColorType[RadiationFreq8]        = 0;
  FieldColorType[RadiationFreq9]        = 0;

  FieldColorType[EmissivityField0]      = 0;
  FieldColorType[EmissivityField1]      = 0;

  FieldColorType[kphHI]                 = 1;
  FieldColorType[PhotoGamma]            = 0;
  FieldColorType[kphHeI]                = 1;
  FieldColorType[kphHeII]               = 1;
  FieldColorType[kdissH2I]              = 1;

  FieldColorType[FieldUndefined]        = 0;

  return SUCCESS;

}
