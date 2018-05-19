/***********************************************************************
/
/  GRID CLASS (ACCESS THE BARYON FIELDS)
/
/  written by: Daniel R. Reynolds
/  date:       October, 2006
/  modified1:
/
/  PURPOSE:
/
/  NOTE:
/
************************************************************************/
 
#include <stdio.h>
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
 
/* function prototypes */
 
int FindField(int f, int farray[], int n);

#ifdef RAD_HYDRO 
 
// Density field 
float* grid::AccessDensity() {
  int DensNum = -1;
  if ((DensNum = FindField(Density, FieldType, NumberOfBaryonFields))<0)
    return NULL;
  return BaryonField[DensNum];
}

// Total Energy field 
float* grid::AccessTotalEnergy() {
  int EnergyNum = -1;
  if ((EnergyNum = FindField(TotalEnergy,FieldType,NumberOfBaryonFields))<0)
    return NULL;
  return BaryonField[EnergyNum];
}

// Gas Energy field 
float* grid::AccessGasEnergy() {
  int GENum = -1;
  if (((GENum = FindField(InternalEnergy,FieldType,NumberOfBaryonFields))<0) 
      || (DualEnergyFormalism == FALSE)) 
    return NULL;
  return BaryonField[GENum];
}

// Velocity1 field 
float* grid::AccessVelocity1() {
  int Vel1Num = -1;
  if ((Vel1Num = FindField(Velocity1, FieldType, NumberOfBaryonFields))<0)
    return NULL;
  return BaryonField[Vel1Num];
}

// Velocity2 field 
float* grid::AccessVelocity2() {
  int Vel2Num = -1;
  if ((Vel2Num = FindField(Velocity2, FieldType, NumberOfBaryonFields))<0) 
    return NULL;
  return BaryonField[Vel2Num];
}

// Velocity3 field 
float* grid::AccessVelocity3() {
  int Vel3Num = -1;
  if ((Vel3Num = FindField(Velocity3, FieldType, NumberOfBaryonFields))<0) 
    return NULL;
  return BaryonField[Vel3Num];
}

// Electron Density field 
float* grid::AccessElectronDensity() {
  int ENum = -1;
  if ((ENum = FindField(ElectronDensity,FieldType,NumberOfBaryonFields))<0) 
    return NULL;
  return BaryonField[ENum];
}

// Hydrogen-I Density field 
float* grid::AccessHIDensity() {
  int HINum = -1;
  if ((HINum = FindField(HIDensity, FieldType, NumberOfBaryonFields))<0) 
    return NULL;
  return BaryonField[HINum];
}

// Hydrogen-II Density field 
float* grid::AccessHIIDensity() {
  int HIINum = -1;
  if ((HIINum = FindField(HIIDensity, FieldType, NumberOfBaryonFields))<0) 
    return NULL;
  return BaryonField[HIINum];
}

// Helium-I Density field 
float* grid::AccessHeIDensity() {
  int HeINum = -1;
  if ((HeINum = FindField(HeIDensity, FieldType, NumberOfBaryonFields))<0) 
    return NULL;
  return BaryonField[HeINum];
}

// Helium-II Density field 
float* grid::AccessHeIIDensity() {
  int HeIINum = -1;
  if ((HeIINum = FindField(HeIIDensity, FieldType, NumberOfBaryonFields))<0) 
    return NULL;
  return BaryonField[HeIINum];
}

// Helium-III Density field 
float* grid::AccessHeIIIDensity() {
  int HeIIINum = -1;
  if ((HeIIINum = FindField(HeIIIDensity, FieldType, NumberOfBaryonFields))<0) 
    return NULL;
  return BaryonField[HeIIINum];
}

// kphHI field 
float* grid::AccessKPhHI() {
  int kphHINum = -1;
  if ((kphHINum = FindField(kphHI, FieldType, NumberOfBaryonFields))<0) 
    return NULL;
  return BaryonField[kphHINum];
}

// PhotoGamma field 
float* grid::AccessPhotoGamma() {
  int PhotoGammaNum = -1;
  if ((PhotoGammaNum = FindField(PhotoGamma, FieldType, NumberOfBaryonFields))<0) 
    return NULL;
  return BaryonField[PhotoGammaNum];
}

// kphHeI field 
float* grid::AccessKPhHeI() {
  int kphHeINum = -1;
  if ((kphHeINum = FindField(kphHeI, FieldType, NumberOfBaryonFields))<0) 
    return NULL;
  return BaryonField[kphHeINum];
}

// kphHeII field 
float* grid::AccessKPhHeII() {
  int kphHeIINum = -1;
  if ((kphHeIINum = FindField(kphHeII, FieldType, NumberOfBaryonFields))<0) 
    return NULL;
  return BaryonField[kphHeIINum];
}

// kdissH2I field 
float* grid::AccessKDissH2I() {
  int kdissH2Inum = -1;
  if ((kdissH2Inum = FindField(kdissH2I, FieldType, NumberOfBaryonFields))<0) 
    return NULL;
  return BaryonField[kdissH2Inum];
}

// Radiation Energy (Grey, or 0th group)
float* grid::AccessRadiationFrequency0() {
  int RadNum = -1;
  if ((RadNum = FindField(RadiationFreq0, FieldType, NumberOfBaryonFields))<0) 
    return NULL;
  return BaryonField[RadNum];
}

// Radiation Energy (1st group)
float* grid::AccessRadiationFrequency1() {
  int RadNum = -1;
  if ((RadNum = FindField(RadiationFreq1, FieldType, NumberOfBaryonFields))<0) 
    return NULL;
  return BaryonField[RadNum];
}

// Radiation Energy (2nd group)
float* grid::AccessRadiationFrequency2() {
  int RadNum = -1;
  if ((RadNum = FindField(RadiationFreq2, FieldType, NumberOfBaryonFields))<0) 
    return NULL;
  return BaryonField[RadNum];
}

// Radiation Energy (3rd group)
float* grid::AccessRadiationFrequency3() {
  int RadNum = -1;
  if ((RadNum = FindField(RadiationFreq3, FieldType, NumberOfBaryonFields))<0) 
    return NULL;
  return BaryonField[RadNum];
}

// Radiation Energy (4th group)
float* grid::AccessRadiationFrequency4() {
  int RadNum = -1;
  if ((RadNum = FindField(RadiationFreq4, FieldType, NumberOfBaryonFields))<0) 
    return NULL;
  return BaryonField[RadNum];
}

// Radiation Energy (5th group)
float* grid::AccessRadiationFrequency5() {
  int RadNum = -1;
  if ((RadNum = FindField(RadiationFreq5, FieldType, NumberOfBaryonFields))<0) 
    return NULL;
  return BaryonField[RadNum];
}

// Radiation Energy (6th group)
float* grid::AccessRadiationFrequency6() {
  int RadNum = -1;
  if ((RadNum = FindField(RadiationFreq6, FieldType, NumberOfBaryonFields))<0) 
    return NULL;
  return BaryonField[RadNum];
}

// Radiation Energy (7th group)
float* grid::AccessRadiationFrequency7() {
  int RadNum = -1;
  if ((RadNum = FindField(RadiationFreq7, FieldType, NumberOfBaryonFields))<0) 
    return NULL;
  return BaryonField[RadNum];
}

// Radiation Energy (8th group)
float* grid::AccessRadiationFrequency8() {
  int RadNum = -1;
  if ((RadNum = FindField(RadiationFreq8, FieldType, NumberOfBaryonFields))<0) 
    return NULL;
  return BaryonField[RadNum];
}

// Radiation Energy (9th group)
float* grid::AccessRadiationFrequency9() {
  int RadNum = -1;
  if ((RadNum = FindField(RadiationFreq9, FieldType, NumberOfBaryonFields))<0) 
    return NULL;
  return BaryonField[RadNum];
}

// Emissivity Field (0th group)
float* grid::AccessEmissivityField0() {
  int EtaNum = -1;
  if ((EtaNum = FindField(EmissivityField0, FieldType, NumberOfBaryonFields))<0) 
    return NULL;
  return BaryonField[EtaNum];
}

// Emissivity Field (1st group)
float* grid::AccessEmissivityField1() {
  int EtaNum = -1;
  if ((EtaNum = FindField(EmissivityField1, FieldType, NumberOfBaryonFields))<0) 
    return NULL;
  return BaryonField[EtaNum];
}

#endif
