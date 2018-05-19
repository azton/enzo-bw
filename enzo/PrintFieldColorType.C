#include <stdio.h>
#include <stdlib.h>

#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"

int PrintFieldColorType(void)
{

  if (debug) {
  printf("Density               %"ISYM"  %"ISYM"\n",      Density              , FieldTypeIsColor(Density)               );
  printf("TotalEnergy           %"ISYM"  %"ISYM"\n",      TotalEnergy          , FieldTypeIsColor(TotalEnergy)           );
  printf("InternalEnergy        %"ISYM"  %"ISYM"\n",      InternalEnergy       , FieldTypeIsColor(InternalEnergy)        );
  printf("Pressure              %"ISYM"  %"ISYM"\n",      Pressure             , FieldTypeIsColor(Pressure)              );
  printf("Velocity1             %"ISYM"  %"ISYM"\n",      Velocity1            , FieldTypeIsColor(Velocity1)             );
  printf("Velocity2             %"ISYM"  %"ISYM"\n",      Velocity2            , FieldTypeIsColor(Velocity2)             );
  printf("Velocity2             %"ISYM"  %"ISYM"\n",      Velocity2            , FieldTypeIsColor(Velocity3)             );
  printf("ElectronDensity       %"ISYM"  %"ISYM"\n",      ElectronDensity      , FieldTypeIsColor(ElectronDensity)       );
  printf("HIDensity             %"ISYM"  %"ISYM"\n",      HIDensity            , FieldTypeIsColor(HIDensity)             );
  printf("HIIDensity            %"ISYM"  %"ISYM"\n",      HIIDensity           , FieldTypeIsColor(HIIDensity)            );
  printf("HeIDensity            %"ISYM"  %"ISYM"\n",      HeIDensity           , FieldTypeIsColor(HeIDensity)            );
  printf("HeIIDensity           %"ISYM"  %"ISYM"\n",      HeIIDensity          , FieldTypeIsColor(HeIIDensity)           );
  printf("HeIIIDensity          %"ISYM"  %"ISYM"\n",      HeIIIDensity         , FieldTypeIsColor(HeIIIDensity)          );
  printf("HMDensity             %"ISYM"  %"ISYM"\n",      HMDensity            , FieldTypeIsColor(HMDensity)             );
  printf("H2IDensity            %"ISYM"  %"ISYM"\n",      H2IDensity           , FieldTypeIsColor(H2IDensity)            );
  printf("H2IIDensity           %"ISYM"  %"ISYM"\n",      H2IIDensity          , FieldTypeIsColor(H2IIDensity)           );
  printf("DIDensity             %"ISYM"  %"ISYM"\n",      DIDensity            , FieldTypeIsColor(DIDensity)             );
  printf("DIIDensity            %"ISYM"  %"ISYM"\n",      DIIDensity           , FieldTypeIsColor(DIIDensity)            );
  printf("HDIDensity            %"ISYM"  %"ISYM"\n",      HDIDensity           , FieldTypeIsColor(HDIDensity)            );
  printf("Metallicity           %"ISYM"  %"ISYM"\n",      Metallicity          , FieldTypeIsColor(Metallicity)           );
  printf("ExtraType0            %"ISYM"  %"ISYM"\n",      ExtraType0           , FieldTypeIsColor(ExtraType0)            );
  printf("ExtraType1            %"ISYM"  %"ISYM"\n",      ExtraType1           , FieldTypeIsColor(ExtraType1)            );
  printf("GravPotential         %"ISYM"  %"ISYM"\n",      GravPotential        , FieldTypeIsColor(GravPotential)         );
  printf("Acceleration0         %"ISYM"  %"ISYM"\n",      Acceleration0        , FieldTypeIsColor(Acceleration0)         );
  printf("Acceleration1         %"ISYM"  %"ISYM"\n",      Acceleration1        , FieldTypeIsColor(Acceleration1)         );
  printf("Acceleration2         %"ISYM"  %"ISYM"\n",      Acceleration2        , FieldTypeIsColor(Acceleration2)         );
  printf("CIDensity             %"ISYM"  %"ISYM"\n",      CIDensity            , FieldTypeIsColor(CIDensity)             );
  printf("CIIDensity            %"ISYM"  %"ISYM"\n",      CIIDensity           , FieldTypeIsColor(CIIDensity)            );
  printf("OIDensity             %"ISYM"  %"ISYM"\n",      OIDensity            , FieldTypeIsColor(OIDensity)             );
  printf("OIIDensity            %"ISYM"  %"ISYM"\n",      OIIDensity           , FieldTypeIsColor(OIIDensity)            );
  printf("SiIDensity            %"ISYM"  %"ISYM"\n",      SiIDensity           , FieldTypeIsColor(SiIDensity)            );
  printf("SiIIDensity           %"ISYM"  %"ISYM"\n",      SiIIDensity          , FieldTypeIsColor(SiIIDensity)           );
  printf("SiIIIDensity          %"ISYM"  %"ISYM"\n",      SiIIIDensity         , FieldTypeIsColor(SiIIIDensity)          );
  printf("CHIDensity            %"ISYM"  %"ISYM"\n",      CHIDensity           , FieldTypeIsColor(CHIDensity)            );
  printf("CH2IDensity           %"ISYM"  %"ISYM"\n",      CH2IDensity          , FieldTypeIsColor(CH2IDensity)           );
  printf("CH3IIDensity          %"ISYM"  %"ISYM"\n",      CH3IIDensity         , FieldTypeIsColor(CH3IIDensity)          );
  printf("C2IDensity            %"ISYM"  %"ISYM"\n",      C2IDensity           , FieldTypeIsColor(C2IDensity)            );
  printf("COIDensity            %"ISYM"  %"ISYM"\n",      COIDensity           , FieldTypeIsColor(COIDensity)            );
  printf("HCOIIDensity          %"ISYM"  %"ISYM"\n",      HCOIIDensity         , FieldTypeIsColor(HCOIIDensity)          );
  printf("OHIDensity            %"ISYM"  %"ISYM"\n",      OHIDensity           , FieldTypeIsColor(OHIDensity)            );
  printf("H2OIDensity           %"ISYM"  %"ISYM"\n",      H2OIDensity          , FieldTypeIsColor(H2OIDensity)           );
  printf("O2IDensity            %"ISYM"  %"ISYM"\n",      O2IDensity           , FieldTypeIsColor(O2IDensity)            );
  printf("Mach                  %"ISYM"  %"ISYM"\n",      Mach                 , FieldTypeIsColor(Mach)                  );
  printf("PreShockTemperature   %"ISYM"  %"ISYM"\n",      PreShockTemperature  , FieldTypeIsColor(PreShockTemperature)   );
  printf("PreShockDensity       %"ISYM"  %"ISYM"\n",      PreShockDensity      , FieldTypeIsColor(PreShockDensity)       );
  printf("CRDensity             %"ISYM"  %"ISYM"\n",      CRDensity            , FieldTypeIsColor(CRDensity)             );
  printf("gParticlePosition     %"ISYM"  %"ISYM"\n",      gParticlePosition    , FieldTypeIsColor(gParticlePosition)     );
  printf("gParticleVelocity     %"ISYM"  %"ISYM"\n",      gParticleVelocity    , FieldTypeIsColor(gParticleVelocity)     );
  printf("gParticleMass         %"ISYM"  %"ISYM"\n",      gParticleMass        , FieldTypeIsColor(gParticleMass)         );
  printf("gParticleAcceleration %"ISYM"  %"ISYM"\n",      gParticleAcceleration, FieldTypeIsColor(gParticleAcceleration) );
  printf("gParticleNumber       %"ISYM"  %"ISYM"\n",      gParticleNumber      , FieldTypeIsColor(gParticleNumber)       );
  printf("gParticleType         %"ISYM"  %"ISYM"\n",      gParticleType        , FieldTypeIsColor(gParticleType)         );
  printf("gParticleAttribute    %"ISYM"  %"ISYM"\n",      gParticleAttribute   , FieldTypeIsColor(gParticleAttribute)    );
  printf("gPotentialField       %"ISYM"  %"ISYM"\n",      gPotentialField      , FieldTypeIsColor(gPotentialField)       );
  printf("gAccelerationField    %"ISYM"  %"ISYM"\n",      gAccelerationField   , FieldTypeIsColor(gAccelerationField)    );
  printf("gGravitatingMassField %"ISYM"  %"ISYM"\n",      gGravitatingMassField, FieldTypeIsColor(gGravitatingMassField) );
  printf("gFlaggingField        %"ISYM"  %"ISYM"\n",      gFlaggingField       , FieldTypeIsColor(gFlaggingField)        );
  printf("gVelocity             %"ISYM"  %"ISYM"\n",      gVelocity            , FieldTypeIsColor(gVelocity)             );
  printf("RadiationFreq0        %"ISYM"  %"ISYM"\n",      RadiationFreq0       , FieldTypeIsColor(RadiationFreq0)        );
  printf("RadiationFreq1        %"ISYM"  %"ISYM"\n",      RadiationFreq1       , FieldTypeIsColor(RadiationFreq1)        );
  printf("RadiationFreq2        %"ISYM"  %"ISYM"\n",      RadiationFreq2       , FieldTypeIsColor(RadiationFreq2)        );
  printf("RadiationFreq3        %"ISYM"  %"ISYM"\n",      RadiationFreq3       , FieldTypeIsColor(RadiationFreq3)        );
  printf("RadiationFreq4        %"ISYM"  %"ISYM"\n",      RadiationFreq4       , FieldTypeIsColor(RadiationFreq4)        );
  printf("RadiationFreq5        %"ISYM"  %"ISYM"\n",      RadiationFreq5       , FieldTypeIsColor(RadiationFreq5)        );
  printf("RadiationFreq6        %"ISYM"  %"ISYM"\n",      RadiationFreq6       , FieldTypeIsColor(RadiationFreq6)        );
  printf("RadiationFreq7        %"ISYM"  %"ISYM"\n",      RadiationFreq7       , FieldTypeIsColor(RadiationFreq7)        );
  printf("RadiationFreq8        %"ISYM"  %"ISYM"\n",      RadiationFreq8       , FieldTypeIsColor(RadiationFreq8)        );
  printf("RadiationFreq9        %"ISYM"  %"ISYM"\n",      RadiationFreq9       , FieldTypeIsColor(RadiationFreq9)        );
  printf("EmissivityField0      %"ISYM"  %"ISYM"\n",      EmissivityField0     , FieldTypeIsColor(EmissivityField0)      );
  printf("EmissivityField1      %"ISYM"  %"ISYM"\n",      EmissivityField1     , FieldTypeIsColor(EmissivityField1)      );
  printf("kphHI                 %"ISYM"  %"ISYM"\n",      kphHI                , FieldTypeIsColor(kphHI)                 );
  printf("PhotoGamma            %"ISYM"  %"ISYM"\n",      PhotoGamma           , FieldTypeIsColor(PhotoGamma)            );
  printf("kphHeI                %"ISYM"  %"ISYM"\n",      kphHeI               , FieldTypeIsColor(kphHeI)                );
  printf("kphHeII               %"ISYM"  %"ISYM"\n",      kphHeII              , FieldTypeIsColor(kphHeII)               );
  printf("kdissH2I              %"ISYM"  %"ISYM"\n",      kdissH2I             , FieldTypeIsColor(kdissH2I)              );
  printf("FieldUndefined        %"ISYM"  %"ISYM"\n",      FieldUndefined       , FieldTypeIsColor(FieldUndefined)        );
  }

  return SUCCESS;

}
