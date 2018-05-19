/*****************************************************************************
 * Copyright 2013 Daniel R. Reynolds                                         *
 *                                                                           *
 * This software is released under the terms of the "Enzo Public License"    *
 * in the accompanying LICENSE file.                                         *
 *****************************************************************************/
/***********************************************************************
   Two-Group (X-ray + UV) Flux-Limited Diffusion Implicit Problem Class
   Parameter output routine

   author: Daniel R. Reynolds
   date:   January 2013

   PURPOSE: Writes internal parameters for problem restart.
************************************************************************/
#ifdef RAD_HYDRO
#include "DualFLD.h"


// function prototypes 
int CosmologyComputeExpansionFactor(FLOAT time, FLOAT *a, FLOAT *dadt);
int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, float *MassUnits, FLOAT Time);


// This routine outputs the restart parameters for the DualFLD solver 
// module.
int DualFLD::WriteParameters(FILE *fptr) {

  // general DualFLD parameters
  fprintf(fptr, "DualFLDChemistry = %"ISYM"\n", Nchem);
  fprintf(fptr, "DualFLDHFraction = %22.16e\n", HFrac);
  fprintf(fptr, "DualFLDMaxDt = %22.16e\n",     maxdt);
  fprintf(fptr, "DualFLDMinDt = %22.16e\n",     mindt);
  fprintf(fptr, "DualFLDInitDt = %22.16e\n",    dtrad);  // set restart initdt to current step
  fprintf(fptr, "DualFLDDtNorm = %22.16e\n",    dtnorm);
  fprintf(fptr, "DualFLDDtGrowth = %22.16e\n",  dtgrowth);
  fprintf(fptr, "DualFLDTheta = %22.16e\n",     theta);
  if (autoScale) {
    fprintf(fptr, "DualFLDAutoScaling = 1\n");
  } else {
    fprintf(fptr, "DualFLDAutoScaling = 0\n");
  }

  // DualFLD X-ray physics parameters
  if (XrayOnly) 
    fprintf(fptr, "DualFLDXrayOnly = 1\n");
  else
    fprintf(fptr, "DualFLDXrayOnly = 0\n");
  if (XrStatic)
    fprintf(fptr, "DualFLDXrayStatic = 1\n");
  else
    fprintf(fptr, "DualFLDXrayStatic = 0\n");
  fprintf(fptr, "DualFLDXraySpectrum = %"ISYM"\n",  XrSpectrum);
  fprintf(fptr, "DualFLDXrayFrequency = %22.16e\n", XrFrequency);
  fprintf(fptr, "DualFLDXrBoundaryX0Faces = %"ISYM" %"ISYM"\n", 
	  XrBdryType[0][0], XrBdryType[0][1]);
  if (rank > 1) {
    fprintf(fptr, "DualFLDXrBoundaryX1Faces = %"ISYM" %"ISYM"\n", 
	    XrBdryType[1][0], XrBdryType[1][1]);
    if (rank > 2) {
      fprintf(fptr, "DualFLDXrBoundaryX2Faces = %"ISYM" %"ISYM"\n", 
	      XrBdryType[2][0], XrBdryType[2][1]);
    }
  }

  // DualFLD X-ray solver parameters
  fprintf(fptr, "DualFLDDtXrayFac = %22.16e\n",        dtfac[1]);
  fprintf(fptr, "DualFLDXrayScaling = %22.16e\n",      XrScale);
  fprintf(fptr, "DualFLDSolToleranceXray = %22.16e\n", sol_tolerance_Xr);
  fprintf(fptr, "DualFLDMaxPCGItersXray = %i\n",       sol_PCGmaxit_Xr);
  fprintf(fptr, "DualFLDMaxMGItersXray = %i\n",        sol_MGmaxit_Xr);
  fprintf(fptr, "DualFLDMGRelaxTypeXray = %i\n",       sol_rlxtype_Xr);    
  fprintf(fptr, "DualFLDMGPreRelaxXray = %i\n",        sol_npre_Xr);
  fprintf(fptr, "DualFLDMGPostRelaxXray = %i\n",       sol_npost_Xr);

  // DualFLD UV physics parameters
  if (UVStatic)
    fprintf(fptr, "DualFLDUVStatic = 1\n");
  else
    fprintf(fptr, "DualFLDUVStatic = 0\n");
  fprintf(fptr, "DualFLDUVSpectrum = %"ISYM"\n",    UVSpectrum);
  fprintf(fptr, "DualFLDUVFrequency = %22.16e\n",   UVFrequency);
  fprintf(fptr, "DualFLDDtUVFac = %22.16e\n",     dtfac[0]);
  fprintf(fptr, "DualFLDUVScaling = %22.16e\n",   UVScale);
  fprintf(fptr, "DualFLDUVBoundaryX0Faces = %"ISYM" %"ISYM"\n", 
	  UVBdryType[0][0], UVBdryType[0][1]);
  if (rank > 1) {
    fprintf(fptr, "DualFLDUVBoundaryX1Faces = %"ISYM" %"ISYM"\n", 
	    UVBdryType[1][0], UVBdryType[1][1]);
    if (rank > 2) {
      fprintf(fptr, "DualFLDUVBoundaryX2Faces = %"ISYM" %"ISYM"\n", 
	      UVBdryType[2][0], UVBdryType[2][1]);
    }
  }

  // DualFLD UV solver parameters
  fprintf(fptr, "DualFLDSolToleranceUV = %22.16e\n",   sol_tolerance_UV);
  fprintf(fptr, "DualFLDMaxMGItersUV = %i\n",          sol_MGmaxit_UV);
  fprintf(fptr, "DualFLDMaxPCGItersUV = %i\n",         sol_PCGmaxit_UV);
  fprintf(fptr, "DualFLDMGRelaxTypeUV = %i\n",         sol_rlxtype_UV);
  fprintf(fptr, "DualFLDMGPreRelaxUV = %i\n",          sol_npre_UV);
  fprintf(fptr, "DualFLDMGPostRelaxUV = %i\n",         sol_npost_UV);

  // if doing an ionization problem (ProblemTypes 210-215),  
  // output additional parameters 
  if ((ProblemType >= 210) && (ProblemType <= 215)) {
    fprintf(fptr, "DualFLDNGammaDotUV = %22.16e\n",   NGammaDotUV);
    fprintf(fptr, "DualFLDNGammaDotXray = %22.16e\n", NGammaDotXr);
    fprintf(fptr, "DualFLDEtaRadius = %22.16e\n",     EtaRadius);
    fprintf(fptr, "DualFLDEtaCenter = %22.16e %22.16e %22.16e\n",
	    EtaCenter[0], EtaCenter[1], EtaCenter[2]);
  }

  // output relevant units: although these aren't required for restart, 
  // cosmology runs never output the units (why?), making data analysis tough
  float TempUnits, RadUnits, MassUnits;
  DenUnits = LenUnits = TempUnits = TimeUnits = VelUnits = MassUnits = 1.0;
  if (GetUnits(&DenUnits, &LenUnits, &TempUnits, &TimeUnits, 
	       &VelUnits, &MassUnits, tnew) != SUCCESS) {
    fprintf(stderr,"DualFLD Evolve: Error in GetUnits.\n");
    return FAIL;
  }
  fprintf(fptr, "DensityUnits = %22.16e\n", DenUnits);
  fprintf(fptr, "LengthUnits = %22.16e\n",  LenUnits);
  fprintf(fptr, "TimeUnits = %22.16e\n",    TimeUnits);

  return SUCCESS;
}

#endif
