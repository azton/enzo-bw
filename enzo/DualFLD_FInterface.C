/*****************************************************************************
 * Copyright 2013 Daniel R. Reynolds                                         *
 *                                                                           *
 * This software is released under the terms of the "Enzo Public License"    *
 * in the accompanying LICENSE file.                                         *
 *****************************************************************************/
/***********************************************************************
   Two-Group (X-ray + UV) Flux-Limited Diffusion Implicit Problem Class
   Fortran interfaces

   author: Daniel R. Reynolds
   date:   January 2013

   PURPOSE: provides C++ interfaces to the relevant Fortran 
            computational kernels
************************************************************************/
#ifdef RAD_HYDRO
#include "DualFLD.h"

// Fortran function prototypes
extern "C" void FORTRAN_NAME(dualfld_setupsystem)(
   Eflt64 *mat, Eflt64 *rhs, float *rhsnorm, float *E, float *kappa, 
   float *eta, float *dt, FLOAT *a, FLOAT *a0, FLOAT *adot, 
   FLOAT *adot0, int *ESpectrum, float *theta, float *aUnits, 
   float *LenUnits, float *LenUnits0, float *EUnits, float *EUnits0, 
   float *NiUnits, float *NiUnits0, float *TUnits, int *rank, float *dx, 
   int *BCTypes, int *N, int *NG, int *faces, int *XrUv, int *ier);

extern "C" void FORTRAN_NAME(dualfld_radiationsource)(
   float *srcUV, float *srcXr, int *ProbType, int *XrOnly, int *UVSpec, 
   float *UVFreq, int *XrSpec, float *XrFreq, float *NGammaDotUV, 
   float *NGammaDotXr, float *EtaRadius, float *EtaCenter, 
   float *LenUnits, int *N, int *NG, float *xB, int *ier);


// C++ interface wrappers
int DualFLD::SetupSystem(float *rhsnorm, int UV_Xr) {
  int ier, ESpectrum;
  float ErUnits, ErUnits0;
  int faces[] = {(OnBdry[0][0]) ? 1 : 0,  (OnBdry[0][1]) ? 1 : 0,
		 (OnBdry[1][0]) ? 1 : 0,  (OnBdry[1][1]) ? 1 : 0,
		 (OnBdry[2][0]) ? 1 : 0,  (OnBdry[2][1]) ? 1 : 0};
  float *E     = U0->GetData(UV_Xr);
  float *kappa = opacity->GetData(UV_Xr);
  float *eta   = extsrc->GetData(UV_Xr);
  if (UV_Xr == 0) {  // Xray
    ESpectrum = XrSpectrum;
    ErUnits   = XrUnits;
    ErUnits0  = XrUnits0;
    FORTRAN_NAME(dualfld_setupsystem)
      (matentries, rhsentries, rhsnorm, E, kappa, eta, &dt, &a, &a0, 
       &adot, &adot0, &ESpectrum, &theta, &aUnits, &LenUnits, &LenUnits0, 
       &ErUnits, &ErUnits0, &NiUnits, &NiUnits0, &TimeUnits, &rank, dx, 
       XrBdryType[0], LocDims, GhDims[0], faces, &UV_Xr, &ier);
  } else {           // UV
    ESpectrum = UVSpectrum;
    ErUnits   = UVUnits;
    ErUnits0  = UVUnits0;
    FORTRAN_NAME(dualfld_setupsystem)
      (matentries, rhsentries, rhsnorm, E, kappa, eta, &dt, &a, &a0, 
       &adot, &adot0, &ESpectrum, &theta, &aUnits, &LenUnits, &LenUnits0, 
       &ErUnits, &ErUnits0, &NiUnits, &NiUnits0, &TimeUnits, &rank, dx, 
       UVBdryType[0], LocDims, GhDims[0], faces, &UV_Xr, &ier);
  }

  // combine the processor-local rhsnorm values together before returning 
  // (since I perform no MPI in F90 modules)
  float rhssum=0.0;
#ifdef USE_MPI
  MPI_Datatype DataType = (sizeof(float) == 4) ? MPI_FLOAT : MPI_DOUBLE;
  MPI_Allreduce(rhsnorm, &rhssum, 1, DataType, MPI_SUM, MPI_COMM_WORLD);
#else
  rhssum = *rhsnorm;
#endif
  *rhsnorm = sqrt(rhssum);

  return(ier);
}


int DualFLD::RadiationSource() {
  int ier, XrOnly=1;
  float lUn    = (LenUnits + LenUnits0)*0.5;
  float *srcXr=NULL, *srcUV=NULL;
  srcXr = extsrc->GetData(0);
  if (!XrayOnly) {
    XrOnly = 0;
    srcUV = extsrc->GetData(1);
  }
  FORTRAN_NAME(dualfld_radiationsource)
    (srcUV, srcXr, &ProblemType, &XrOnly, &UVSpectrum, &UVFrequency, 
     &XrSpectrum, &XrFrequency, &NGammaDotUV, &NGammaDotXr, &EtaRadius, 
     EtaCenter, &lUn, LocDims, GhDims[0], EdgeVals[0], &ier);
  return(ier);
}

#endif
