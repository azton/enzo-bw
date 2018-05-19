/*****************************************************************************
 * Copyright 2013 Daniel R. Reynolds                                         *
 *                                                                           *
 * This software is released under the terms of the "Enzo Public License"    *
 * in the accompanying LICENSE file.                                         *
 *****************************************************************************/
/***********************************************************************
   Two-Group (X-ray + UV) Flux-Limited Diffusion Implicit Problem Class
   Constructor routine

   author: Daniel R. Reynolds
   date:   January 2013

   PURPOSE: Initializes all values to illegal numbers, and sets all 
            arrays to NULL;  Requires call to Initialize to actually 
            set up these values.
************************************************************************/
#ifdef RAD_HYDRO
#include "DualFLD.h"


DualFLD::DualFLD() {

  int dim, face, i;

  // initialize total RT time to zero
  RTtime    = 0.0;
  HYPREtime = 0.0;

  // initialize HYPRE values to -1/NULL
  mattype          = -1;
  stSize           = -1;
  grid             = NULL;
  stencil          = NULL;
  sol_tolerance_Xr = -1.0;
  sol_tolerance_UV = -1.0;
  sol_MGmaxit_Xr   = -1;
  sol_MGmaxit_UV   = -1;
  sol_PCGmaxit_Xr  = -1;
  sol_PCGmaxit_UV  = -1;
  sol_rlxtype_Xr   = -1;
  sol_rlxtype_UV   = -1;
  sol_npre_Xr      = -1;
  sol_npre_UV      = -1;
  sol_npost_Xr     = -1;
  sol_npost_UV     = -1;
  sol_printl       = -1;
  sol_log          = -1;
  totIters         = -1;
  for (dim=0; dim<3; dim++) {
    for (face=0; face<2; face++)
      SolvIndices[dim][face] = 0;
    SolvOff[dim] = 0;
  }

  // initialize problem grid information to -1/NULL
  rank = -1;
  for (dim=0; dim<3; dim++) {
    layout[dim]   = 1;  // initialize for single-processor run
    location[dim] = 0;  // initialize for single-processor run
    LocDims[dim]  = 1;
    ArrDims[dim]  = 1;
    dx[dim]       = 1.0;
    for (face=0; face<2; face++) {
      OnBdry[dim][face]   = false;
      NBors[dim][face]    = MPI_PROC_NULL;
      GhDims[dim][face]   = 0;
      XrBdryType[dim][face] = -1;
      UVBdryType[dim][face] = -1;
      EdgeVals[dim][face] = -1.0;
      XrBdryVals[dim][face] = NULL;
      UVBdryVals[dim][face] = NULL;
    }
  }
  
  // initialize time-stepping related data to -1/NULL
  maxdt    = 1.0e20;
  mindt    = 0.0;
  initdt   = 1.0e20;
  dtfac[0] = 1.0e20;
  dtfac[1] = 1.0e20;
  dtnorm   = 0.0;
  dtgrowth = 1.1;
  tnew     = -1.0;
  told     = -1.0;
  dt       = -1.0;
  dtrad    = -1.0;
  theta    = -1.0;
  sol      = NULL;
  U0       = NULL;
  extsrc   = NULL;
  opacity  = NULL;
  
  // initialize problem defining data 
  XrayOnly    = false;
  XrStatic    = false;
  UVStatic    = false;
  a           = 1.0;
  a0          = 1.0;
  adot        = 0.0;
  adot0       = 0.0;
  aUnits      = 1.0;
  UVScale     = 1.0;
  XrScale     = 1.0;
  autoScale   = true;
  UVUnits     = 1.0;
  XrUnits0    = 1.0;
  NiUnits     = 1.0;
  NiUnits0    = 1.0;
  DenUnits    = 1.0;
  DenUnits0   = 1.0;
  LenUnits    = 1.0;
  LenUnits0   = 1.0;
  TimeUnits   = 1.0;
  VelUnits    = 1.0;
  Nchem       = -1;
  HFrac       = -1.0;
  UVSpectrum  = -1;
  UVFrequency = 0.0;
  XrSpectrum  = -1;
  XrFrequency = 0.0;
  NGammaDotUV = 0.0;
  NGammaDotXr = 0.0;

  for (i=0; i<7; i++)  RadIntUV[i] = 0.0;
  for (i=0; i<7; i++)  RadIntXr[i] = 0.0;

  // initialize linear solver/Jacobian arrays to NULL
  matentries = NULL;
  rhsentries = NULL;
  HYPREbuff  = NULL;

  // initialize HYPRE structures to NULL
  P      = NULL;
  rhsvec = NULL;
  solvec = NULL;
}

#endif
