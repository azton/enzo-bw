#include <stdio.h>
#include <stdlib.h>

#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"

int PrintFieldDensityType(void)
{

  if (debug) {
  printf("Density               %"ISYM"  %"ISYM"\n",      Density              , FieldTypeIsDensity(Density)               );
  printf("TotalEnergy           %"ISYM"  %"ISYM"\n",      TotalEnergy          , FieldTypeIsDensity(TotalEnergy)           );
  printf("InternalEnergy        %"ISYM"  %"ISYM"\n",      InternalEnergy       , FieldTypeIsDensity(InternalEnergy)        );
  printf("Pressure              %"ISYM"  %"ISYM"\n",      Pressure             , FieldTypeIsDensity(Pressure)              );
  printf("Velocity1             %"ISYM"  %"ISYM"\n",      Velocity1            , FieldTypeIsDensity(Velocity1)             );
  printf("Velocity2             %"ISYM"  %"ISYM"\n",      Velocity2            , FieldTypeIsDensity(Velocity2)             );
  printf("Velocity2             %"ISYM"  %"ISYM"\n",      Velocity2            , FieldTypeIsDensity(Velocity3)             );
  printf("ElectronDensity       %"ISYM"  %"ISYM"\n",      ElectronDensity      , FieldTypeIsDensity(ElectronDensity)       );
  printf("HIDensity             %"ISYM"  %"ISYM"\n",      HIDensity            , FieldTypeIsDensity(HIDensity)             );
  printf("HIIDensity            %"ISYM"  %"ISYM"\n",      HIIDensity           , FieldTypeIsDensity(HIIDensity)            );
  printf("HeIDensity            %"ISYM"  %"ISYM"\n",      HeIDensity           , FieldTypeIsDensity(HeIDensity)            );
  printf("HeIIDensity           %"ISYM"  %"ISYM"\n",      HeIIDensity          , FieldTypeIsDensity(HeIIDensity)           );
  printf("HeIIIDensity          %"ISYM"  %"ISYM"\n",      HeIIIDensity         , FieldTypeIsDensity(HeIIIDensity)          );
  printf("HMDensity             %"ISYM"  %"ISYM"\n",      HMDensity            , FieldTypeIsDensity(HMDensity)             );
  printf("H2IDensity            %"ISYM"  %"ISYM"\n",      H2IDensity           , FieldTypeIsDensity(H2IDensity)            );
  printf("H2IIDensity           %"ISYM"  %"ISYM"\n",      H2IIDensity          , FieldTypeIsDensity(H2IIDensity)           );
  printf("DIDensity             %"ISYM"  %"ISYM"\n",      DIDensity            , FieldTypeIsDensity(DIDensity)             );
  printf("DIIDensity            %"ISYM"  %"ISYM"\n",      DIIDensity           , FieldTypeIsDensity(DIIDensity)            );
  printf("HDIDensity            %"ISYM"  %"ISYM"\n",      HDIDensity           , FieldTypeIsDensity(HDIDensity)            );
  printf("Metallicity           %"ISYM"  %"ISYM"\n",      Metallicity          , FieldTypeIsDensity(Metallicity)           );
  printf("ExtraType0            %"ISYM"  %"ISYM"\n",      ExtraType0           , FieldTypeIsDensity(ExtraType0)            );
  printf("ExtraType1            %"ISYM"  %"ISYM"\n",      ExtraType1           , FieldTypeIsDensity(ExtraType1)            );
  printf("GravPotential         %"ISYM"  %"ISYM"\n",      GravPotential        , FieldTypeIsDensity(GravPotential)         );
  printf("Acceleration0         %"ISYM"  %"ISYM"\n",      Acceleration0        , FieldTypeIsDensity(Acceleration0)         );
  printf("Acceleration1         %"ISYM"  %"ISYM"\n",      Acceleration1        , FieldTypeIsDensity(Acceleration1)         );
  printf("Acceleration2         %"ISYM"  %"ISYM"\n",      Acceleration2        , FieldTypeIsDensity(Acceleration2)         );
  printf("CIDensity             %"ISYM"  %"ISYM"\n",      CIDensity            , FieldTypeIsDensity(CIDensity)             );
  printf("CIIDensity            %"ISYM"  %"ISYM"\n",      CIIDensity           , FieldTypeIsDensity(CIIDensity)            );
  printf("OIDensity             %"ISYM"  %"ISYM"\n",      OIDensity            , FieldTypeIsDensity(OIDensity)             );
  printf("OIIDensity            %"ISYM"  %"ISYM"\n",      OIIDensity           , FieldTypeIsDensity(OIIDensity)            );
  printf("SiIDensity            %"ISYM"  %"ISYM"\n",      SiIDensity           , FieldTypeIsDensity(SiIDensity)            );
  printf("SiIIDensity           %"ISYM"  %"ISYM"\n",      SiIIDensity          , FieldTypeIsDensity(SiIIDensity)           );
  printf("SiIIIDensity          %"ISYM"  %"ISYM"\n",      SiIIIDensity         , FieldTypeIsDensity(SiIIIDensity)          );
  printf("CHIDensity            %"ISYM"  %"ISYM"\n",      CHIDensity           , FieldTypeIsDensity(CHIDensity)            );
  printf("CH2IDensity           %"ISYM"  %"ISYM"\n",      CH2IDensity          , FieldTypeIsDensity(CH2IDensity)           );
  printf("CH3IIDensity          %"ISYM"  %"ISYM"\n",      CH3IIDensity         , FieldTypeIsDensity(CH3IIDensity)          );
  printf("C2IDensity            %"ISYM"  %"ISYM"\n",      C2IDensity           , FieldTypeIsDensity(C2IDensity)            );
  printf("COIDensity            %"ISYM"  %"ISYM"\n",      COIDensity           , FieldTypeIsDensity(COIDensity)            );
  printf("HCOIIDensity          %"ISYM"  %"ISYM"\n",      HCOIIDensity         , FieldTypeIsDensity(HCOIIDensity)          );
  printf("OHIDensity            %"ISYM"  %"ISYM"\n",      OHIDensity           , FieldTypeIsDensity(OHIDensity)            );
  printf("H2OIDensity           %"ISYM"  %"ISYM"\n",      H2OIDensity          , FieldTypeIsDensity(H2OIDensity)           );
  printf("O2IDensity            %"ISYM"  %"ISYM"\n",      O2IDensity           , FieldTypeIsDensity(O2IDensity)            );
  printf("Mach                  %"ISYM"  %"ISYM"\n",      Mach                 , FieldTypeIsDensity(Mach)                  );
  printf("PreShockTemperature   %"ISYM"  %"ISYM"\n",      PreShockTemperature  , FieldTypeIsDensity(PreShockTemperature)   );
  printf("PreShockDensity       %"ISYM"  %"ISYM"\n",      PreShockDensity      , FieldTypeIsDensity(PreShockDensity)       );
  printf("CRDensity             %"ISYM"  %"ISYM"\n",      CRDensity            , FieldTypeIsDensity(CRDensity)             );
  printf("gParticlePosition     %"ISYM"  %"ISYM"\n",      gParticlePosition    , FieldTypeIsDensity(gParticlePosition)     );
  printf("gParticleVelocity     %"ISYM"  %"ISYM"\n",      gParticleVelocity    , FieldTypeIsDensity(gParticleVelocity)     );
  printf("gParticleMass         %"ISYM"  %"ISYM"\n",      gParticleMass        , FieldTypeIsDensity(gParticleMass)         );
  printf("gParticleAcceleration %"ISYM"  %"ISYM"\n",      gParticleAcceleration, FieldTypeIsDensity(gParticleAcceleration) );
  printf("gParticleNumber       %"ISYM"  %"ISYM"\n",      gParticleNumber      , FieldTypeIsDensity(gParticleNumber)       );
  printf("gParticleType         %"ISYM"  %"ISYM"\n",      gParticleType        , FieldTypeIsDensity(gParticleType)         );
  printf("gParticleAttribute    %"ISYM"  %"ISYM"\n",      gParticleAttribute   , FieldTypeIsDensity(gParticleAttribute)    );
  printf("gPotentialField       %"ISYM"  %"ISYM"\n",      gPotentialField      , FieldTypeIsDensity(gPotentialField)       );
  printf("gAccelerationField    %"ISYM"  %"ISYM"\n",      gAccelerationField   , FieldTypeIsDensity(gAccelerationField)    );
  printf("gGravitatingMassField %"ISYM"  %"ISYM"\n",      gGravitatingMassField, FieldTypeIsDensity(gGravitatingMassField) );
  printf("gFlaggingField        %"ISYM"  %"ISYM"\n",      gFlaggingField       , FieldTypeIsDensity(gFlaggingField)        );
  printf("gVelocity             %"ISYM"  %"ISYM"\n",      gVelocity            , FieldTypeIsDensity(gVelocity)             );
  printf("RadiationFreq0        %"ISYM"  %"ISYM"\n",      RadiationFreq0       , FieldTypeIsDensity(RadiationFreq0)        );
  printf("RadiationFreq1        %"ISYM"  %"ISYM"\n",      RadiationFreq1       , FieldTypeIsDensity(RadiationFreq1)        );
  printf("RadiationFreq2        %"ISYM"  %"ISYM"\n",      RadiationFreq2       , FieldTypeIsDensity(RadiationFreq2)        );
  printf("RadiationFreq3        %"ISYM"  %"ISYM"\n",      RadiationFreq3       , FieldTypeIsDensity(RadiationFreq3)        );
  printf("RadiationFreq4        %"ISYM"  %"ISYM"\n",      RadiationFreq4       , FieldTypeIsDensity(RadiationFreq4)        );
  printf("RadiationFreq5        %"ISYM"  %"ISYM"\n",      RadiationFreq5       , FieldTypeIsDensity(RadiationFreq5)        );
  printf("RadiationFreq6        %"ISYM"  %"ISYM"\n",      RadiationFreq6       , FieldTypeIsDensity(RadiationFreq6)        );
  printf("RadiationFreq7        %"ISYM"  %"ISYM"\n",      RadiationFreq7       , FieldTypeIsDensity(RadiationFreq7)        );
  printf("RadiationFreq8        %"ISYM"  %"ISYM"\n",      RadiationFreq8       , FieldTypeIsDensity(RadiationFreq8)        );
  printf("RadiationFreq9        %"ISYM"  %"ISYM"\n",      RadiationFreq9       , FieldTypeIsDensity(RadiationFreq9)        );
  printf("EmissivityField0      %"ISYM"  %"ISYM"\n",      EmissivityField0     , FieldTypeIsDensity(EmissivityField0)      );
  printf("EmissivityField1      %"ISYM"  %"ISYM"\n",      EmissivityField1     , FieldTypeIsDensity(EmissivityField1)      );
  printf("kphHI                 %"ISYM"  %"ISYM"\n",      kphHI                , FieldTypeIsDensity(kphHI)                 );
  printf("PhotoGamma            %"ISYM"  %"ISYM"\n",      PhotoGamma           , FieldTypeIsDensity(PhotoGamma)            );
  printf("kphHeI                %"ISYM"  %"ISYM"\n",      kphHeI               , FieldTypeIsDensity(kphHeI)                );
  printf("kphHeII               %"ISYM"  %"ISYM"\n",      kphHeII              , FieldTypeIsDensity(kphHeII)               );
  printf("kdissH2I              %"ISYM"  %"ISYM"\n",      kdissH2I             , FieldTypeIsDensity(kdissH2I)              );
  printf("FieldUndefined        %"ISYM"  %"ISYM"\n",      FieldUndefined       , FieldTypeIsDensity(FieldUndefined)        );
  }

  return SUCCESS;

}
