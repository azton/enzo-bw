/*****************************************************************************
 * Copyright 2013 Daniel R. Reynolds                                         *
 *                                                                           *
 * This software is released under the terms of the "Enzo Public License"    *
 * in the accompanying LICENSE file.                                         *
 *****************************************************************************/
/***********************************************************************
   Two-Group (X-ray + UV) Flux-Limited Diffusion Implicit Problem Class
   Evolve Routine for split radiation-chemistry system

   author: Daniel R. Reynolds
   date:   January 2013

   PURPOSE: Computes the UV and X-ray opacities over the domain 

   NOTE: These are normalized by NiUnits, and must be converted to 
         the appropriate redshift-dependent values when used, e.g.
              opacity_true = opacity*NiUnits.

   RETURNS: SUCCESS or FAIL
************************************************************************/
#ifdef RAD_HYDRO

#include "DualFLD.h"

int DualFLD::ComputeOpacity(float *HI, float *HeI, float *HeII) {

  // local variables
  float *kappaXr=NULL, *kappaUV=NULL;
  kappaXr = opacity->GetData(0);
  if (!XrayOnly)  kappaUV = opacity->GetData(1);
  float cHI_Xr, cHeI_Xr, cHeII_Xr, cHI_UV, cHeI_UV, cHeII_UV;
  int i, size=1;
  for (i=0; i<rank; i++)  size *= ArrDims[i];

  // initialize outputs to zero
  opacity->constant(0.0);

  // set shortcuts
  cHI_Xr   = RadIntXr[1] / RadIntXr[0];
  cHeI_Xr  = RadIntXr[3] / RadIntXr[0] / 4.0;
  cHeII_Xr = RadIntXr[5] / RadIntXr[0] / 4.0;
  cHI_UV   = RadIntUV[1] / RadIntUV[0];
  cHeI_UV  = RadIntUV[3] / RadIntUV[0] / 4.0;
  cHeII_UV = RadIntUV[5] / RadIntUV[0] / 4.0;

  // compute opacities depending on chemistry
  if (Nchem == 1) {           // hydrogen only

    #pragma omp parallel default(shared)
    {
      #pragma omp for schedule(static)
      for (i=0; i<size; i++)  kappaXr[i] = HI[i]*cHI_Xr;
      if (!XrayOnly)
        #pragma omp for schedule(static)
	for (i=0; i<size; i++)  kappaUV[i] = HI[i]*cHI_UV;
    }

  } else if (Nchem == 3) {    // hydrogen + helium

    #pragma omp parallel default(shared)
    {
      #pragma omp for schedule(static)
      for (i=0; i<size; i++) 
	kappaXr[i] = HI[i]*cHI_Xr + HeI[i]*cHeI_Xr + HeII[i]*cHeII_Xr;
      if (!XrayOnly)
        #pragma omp for schedule(static)
	for (i=0; i<size; i++) 
	  kappaUV[i] = HI[i]*cHI_UV + HeI[i]*cHeI_UV + HeII[i]*cHeII_UV;
    }

  } else {                    // illegal Nchem
    fprintf(stderr,"DualFLD_ComputeOpacity ERROR: illegal Nchem = %"ISYM
	    ", requires 1 or 3.\n", Nchem);
    return FAIL;

  }

  return SUCCESS;
}

#endif
