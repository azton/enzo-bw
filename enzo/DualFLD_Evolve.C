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

   PURPOSE: Evolves the FLD system forward in time to "catch up" with 
            the hydrodynamics time.

   RETURNS: SUCCESS or FAIL
************************************************************************/
#ifdef RAD_HYDRO

#define SOLUTION_DIAGNOSTICS

#ifdef USE_PAT
#include <pat_api.h>
#endif /* USE_PAT */
 
#include "DualFLD.h"
#include "CosmologyParameters.h"


// function prototypes 
int CosmologyComputeExpansionFactor(FLOAT time, FLOAT *a, FLOAT *dadt);
int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, float *MassUnits, FLOAT Time);
int RadiationGetUnits(float *RadiationUnits, FLOAT Time);


// This routine evolves the radiation equations in an operator-split fashion:
// we take radiation time steps to catch up with the hydrodynamic time, i.e. 
// (dt_rad <= dt_hydro). Both the UV and X-ray radiation are computed with 
// the same step size.  Since chemistry and gas cooling are handled externally 
// to this module, we generally prefer (dt_rad == dt_hydro), but subcycling is
// allowed for robustness when handling fast transients (e.g. the first step
// after a source turns on).  As a result, prior to completion this routine
// updates the maximum time step the overall Grid module can take so that
// subsequent steps are taken with the same sizes.
int DualFLD::Evolve(HierarchyEntry *ThisGrid, float dthydro) {

  // Only continue if we own this grid
  if (MyProcessorNumber != ThisGrid->GridData->ReturnProcessorNumber())
    return SUCCESS;

  // in case MPI is not included
#ifndef MPI_INT
  int MPI_COMM_WORLD = 0;
#endif

  // start MPI timer for overall solver
#ifdef USE_MPI
  float stime = MPI_Wtime();
#else
  float stime = 0.0;
#endif

  ////////////////////////////////////
  // Problem Setup Phase

  if (debug)  printf("\n DualFLD Evolve:\n");

  // Set pointers to each variable (zero out fluid energy correction 
  // since that is not internal to Enzo)
  float *UVRadiation=NULL, *XrRadiation=NULL;
  XrRadiation = ThisGrid->GridData->AccessRadiationFrequency0();
  if (XrRadiation == NULL) {
    fprintf(stderr,"DualFLD Evolve: could not obtain X-ray radiation\n");
    return FAIL; }
  if (!XrayOnly) {
    UVRadiation = ThisGrid->GridData->AccessRadiationFrequency1();
    if (UVRadiation == NULL) {
      fprintf(stderr,"DualFLD Evolve: could not obtain UV radiation\n");
      return FAIL; }
  }

  // Get current time from Grid
  tnew = ThisGrid->GridData->ReturnTime();

  // initialize external sources to 0
  extsrc->constant(0.0);

  // access/fill radiation source array (if provided externally)
  int i, size=1;
  float *XrSrc=NULL, *UVSrc=NULL;
  for (i=0; i<rank; i++)  size *= ArrDims[i];
  XrSrc = extsrc->GetData(0);
  if (XrSrc == NULL) {
    fprintf(stderr,"DualFLD Evolve: could not access X-ray source\n");
    return FAIL; }
  if (!XrayOnly) {
    UVSrc = extsrc->GetData(1);
    if (UVSrc == NULL) {
      fprintf(stderr,"DualFLD Evolve: could not access UV source\n");
      return FAIL; }
  }
  int eta_set = 0;
#ifdef EMISSIVITY
  // if using external Emissivity field source, copy into extsrc
  if (StarMakerEmissivityField > 0) {
    // // access external emissivity fields
    // float *UVSource=NULL, *XrSource=NULL;
    // XrSource = ThisGrid->GridData->AccessEmissivityField0();
    // if (XrSource == NULL) {
    //   fprintf(stderr,"DualFLD Evolve: could not access X-ray emissivity field\n");
    //   return FAIL; }
    // if (!XrayOnly) {
    //   UVSource = ThisGrid->GridData->AccessEmissivityField1();
    //   if (UVSource == NULL) {
    // 	fprintf(stderr,"DualFLD Evolve: could not access UV emissivity field\n");
    // 	return FAIL; }
    // }
    // // copy data
    // for (i=0; i<size; i++)  XrSrc[i] = XrSource[i];
    // if (!XrayOnly)
    //   for (i=0; i<size; i++)  UVSrc[i] = UVSource[i];

    // access external emissivity fields
    float *UVSource=NULL;
    UVSource = ThisGrid->GridData->AccessEmissivityField0();
    // copy data
    float factor = Xr_parameter / UV_parameter;
    for (i=0; i<size; i++)  XrSrc[i] = UVSource[i] * factor;
    if (!XrayOnly)
      for (i=0; i<size; i++)  UVSrc[i] = UVSource[i];

    eta_set = 1;
#ifdef SOLUTION_DIAGNOSTICS
    float Xr_norm=0.0, UV_norm=0.0;
    if (XrayOnly) {
      Xr_norm = extsrc->rmsnorm_component(0);
      if (debug) printf("   emissivity norm:  X-ray = %g\n", Xr_norm);
    } else {
      Xr_norm = extsrc->rmsnorm_component(0);
      UV_norm = extsrc->rmsnorm_component(1);
      if (debug) printf("   emissivity norms:  X-ray = %g,  UV = %g\n",
			Xr_norm, UV_norm);
    }
#endif
  }
#endif

  // attach arrays to U0 vector
  U0->SetData(0, XrRadiation);
  if (!XrayOnly)
    U0->SetData(1, UVRadiation);

  // rescale Enzo units with input scalings to non-dimensionalize within solver
  U0->scale_component(0, 1.0/XrScale);
  if (!XrayOnly)
    U0->scale_component(1, 1.0/UVScale);

  // have U0 begin communication of neighbor information
  if (U0->exchange_start() != SUCCESS) {
    fprintf(stderr,"DualFLD Evolve: vector exchange_start error\n");
    return FAIL; }

  // update internal Enzo units for current times
  float TempUnits, RadUnits, MassUnits;
  DenUnits = LenUnits = TempUnits = TimeUnits = VelUnits = MassUnits = 1.0;
  if (GetUnits(&DenUnits, &LenUnits, &TempUnits, &TimeUnits, 
	       &VelUnits, &MassUnits, told) != SUCCESS) {
    fprintf(stderr,"DualFLD Evolve: Error in GetUnits.\n");
    return FAIL;
  }
  if (RadiationGetUnits(&RadUnits, told) != SUCCESS) {
    fprintf(stderr,"DualFLD Evolve: Error in RadiationGetUnits.\n");
    return FAIL;
  }
  float mp = 1.67262171e-24;   // Mass of a proton [g]
  UVUnits0 = RadUnits*UVScale;
  XrUnits0 = RadUnits*XrScale;
  NiUnits0 = (Nchem == 0) ? 1.0 : DenUnits/mp;

  DenUnits = LenUnits = TempUnits = TimeUnits = VelUnits = MassUnits = 1.0;
  if (GetUnits(&DenUnits, &LenUnits, &TempUnits, &TimeUnits, 
	       &VelUnits, &MassUnits, tnew) != SUCCESS) {
    fprintf(stderr,"DualFLD Evolve: Error in GetUnits.\n");
    return FAIL;
  }
  if (RadiationGetUnits(&RadUnits, tnew) != SUCCESS) {
    fprintf(stderr,"DualFLD Evolve: Error in RadiationGetUnits.\n");
    return FAIL;
  }
  UVUnits = RadUnits*UVScale;
  XrUnits = RadUnits*XrScale;
  NiUnits = (Nchem == 0) ? 1.0 : DenUnits/mp;

#ifdef SOLUTION_DIAGNOSTICS
  float UTypVals[2], UMaxVals[2];
  // output typical/maximum values
  if (XrayOnly) {
    UTypVals[0] = U0->rmsnorm_component(0);
    UMaxVals[0] = U0->infnorm_component(0);

    if (debug) {
      printf("   current internal (physical) quantities:\n");
      printf("    Xray rms = %13.7e (%8.2e), max = %13.7e (%8.2e)\n",
	     UTypVals[0], UTypVals[0]*XrUnits, 
	     UMaxVals[0], UMaxVals[0]*XrUnits);
    }
  } else {
    UTypVals[0] = U0->rmsnorm_component(0);
    UTypVals[1] = U0->rmsnorm_component(1);
    UMaxVals[0] = U0->infnorm_component(0);
    UMaxVals[1] = U0->infnorm_component(1);

    if (debug) {
      printf("   current internal (physical) quantities:\n");
      printf("    Xray rms = %13.7e (%8.2e), max = %13.7e (%8.2e)\n",
	     UTypVals[0], UTypVals[0]*UVUnits, UMaxVals[0], UMaxVals[0]*UVUnits);
      printf("      UV rms = %13.7e (%8.2e), max = %13.7e (%8.2e)\n",
	     UTypVals[1], UTypVals[1]*XrUnits, UMaxVals[1], UMaxVals[1]*XrUnits);
    }
  }
#endif

  // if autoScale enabled, determine scaling factor updates here
  float ScaleCorrTol = 1.e-2;
  float XrScaleCorr = 1.0;
  float UVScaleCorr = 1.0;
  if (autoScale) {

#ifndef SOLUTION_DIAGNOSTICS
    float UTypVals[2], UMaxVals[2];
    if (XrayOnly) {
      UTypVals[0] = U0->rmsnorm_component(0);
      UMaxVals[0] = U0->infnorm_component(0);
    } else {
      UTypVals[0] = U0->rmsnorm_component(0);
      UTypVals[1] = U0->rmsnorm_component(1);
      UMaxVals[0] = U0->infnorm_component(0);
      UMaxVals[1] = U0->infnorm_component(1);
    }
#endif
    if ((UMaxVals[0] - UTypVals[0]) > ScaleCorrTol*UMaxVals[0])
      XrScaleCorr = UMaxVals[0];
    if ((UMaxVals[1] - UTypVals[1]) > ScaleCorrTol*UMaxVals[1])
      UVScaleCorr = UMaxVals[1];
  }

  // have U0 finish communication of neighbor information
  if (U0->exchange_end() != SUCCESS) {
    fprintf(stderr,"DualFLD Evolve: vector exchange_end error\n");
    return FAIL;
  }


  ////////////////////////////////////
  // Problem Solve Phase

  // internal time-stepping loop to catch up with Hydro time
  float stime2, ftime2;   // radiation timer
  int radstep, radstop;   // subcycle iterators
  float end_time = tnew + dthydro;
  radstop = 0;
  for (radstep=0; radstep<=100; radstep++) {
      
    // start MPI timer for radiation solver
    stime2 = MPI_Wtime();
  
    // update time-step information
    told = tnew;

    // keep trying time steps until radiation solver succeeds. 
    // Note: if we reach the minimum time step size, RadStep will FAIL
    int recompute_step = 0;
    do {

      // update time-step information.  Note: dtrad was set on previous 
      // iteration of solver, or by user input for first iteration
      tnew = told + dtrad;
      if ((tnew - end_time)/end_time > -1.0e-14) {   // don't pass synchronization time
	tnew = end_time;
	radstop = 1;
      }
      dt = tnew - told;
      if (debug) {
	printf("\n subcycled rad %"ISYM": dt=%7.1e, t=%7.1e (hydro dt=%7.1e, t=%7.1e)\n",
	       radstep,dt,tnew,dthydro,end_time);
	printf(" ----------------------------------------------------------------------\n");
      }

      // take a Xray radiation step (if not static)
      if (XrStatic) {
	if (debug)  printf("  Xray: field is static, so just copying old values\n");
	sol->copy_component(U0, 0);
      } else {
	if (debug)  printf("  Xray:");
	recompute_step = this->RadStep(ThisGrid, eta_set, 0);
      }

      // take a UV radiation step (if Xray step succeeded, and UV enabled)
      if (!recompute_step && !XrayOnly) {
	if (UVStatic) {
	  if (debug)  printf("  UV: field is static, so just copying old values\n");
	  sol->copy_component(U0, 1);
	} else {
	  if (debug)  printf("  UV:");
	  recompute_step = this->RadStep(ThisGrid, eta_set, 1);
	}
      }
      
      if (debug)  
	printf(" ======================================================================\n\n");

      // if either radiation step was unsuccessful, update dtrad and try again
      if (recompute_step) {
	dtrad = max(dtrad*0.5, mindt);
	radstop = 0;
	tnew = told;
      }
      
    } while(recompute_step);
    
    // stop MPI timer for radiation solver, increment total
    ftime2 = MPI_Wtime();
    HYPREtime += ftime2-stime2;
    
    // update the radiation time step size for next time step
    //   (limit growth at each cycle)
    float dt_est;
    if (XrayOnly)
      dt_est = this->ComputeTimeStep(U0,sol,0);
    else 
      dt_est = this->ComputeTimeStep(U0,sol,2);
    dtrad = min(dt_est, dtgrowth*dtrad);

    // update Enzo radiation fields with new values
    U0->copy(sol);

    // have U0 communicate neighbor information
    if (U0->exchange() != SUCCESS) {
      fprintf(stderr,"DualFLD Evolve: vector exchange error\n");
      return FAIL;
    }

    // break out of time-stepping loop if we've reached the end
    if (radstop)  break;
	
  } // end outer radiation time-stepping loop


  ////////////////////////////////////
  // Problem Cleanup and Preparation for Next Call Phase

  // if chemistry and cooling handled elsewhere, fill the rates
  if (RadiativeCooling) {
    float *rho        = ThisGrid->GridData->AccessDensity();
    float *HI         = ThisGrid->GridData->AccessHIDensity();
    float *HeI        = ThisGrid->GridData->AccessHeIDensity();
    float *HeII       = ThisGrid->GridData->AccessHeIIDensity();
    float *phHI       = ThisGrid->GridData->AccessKPhHI();
    float *phHeI      = ThisGrid->GridData->AccessKPhHeI();
    float *phHeII     = ThisGrid->GridData->AccessKPhHeII();
    float *photogamma = ThisGrid->GridData->AccessPhotoGamma();
    float *dissH2I    = ThisGrid->GridData->AccessKDissH2I();
    this->FillRates(sol, rho, HI, HeI, HeII, phHI, phHeI, 
		    phHeII, photogamma, dissH2I);
  }

  // update the radiation time step size for next time step
  if (dtrad != huge_number)
    ThisGrid->GridData->SetMaxRadiationDt(dtrad);

  // scale back to Enzo units
  U0->scale_component(0,XrScale);
  if (!XrayOnly)  U0->scale_component(1,UVScale);

  // update scaling factors to account for new values
  if (autoScale) {
    XrScale *= XrScaleCorr;
    if (!XrayOnly)  UVScale *= UVScaleCorr;
  }

  // stop MPI timer, add to cumulative clock, output to stdout
#ifdef USE_MPI
  float ftime = MPI_Wtime();
#else
  float ftime = 0.0;
#endif
  RTtime += ftime-stime;
  if (debug)  printf("RadHydro cumulative time = %g (HYPRE = %g)\n\n",
		     RTtime, HYPREtime);

  return SUCCESS;
}



// This routine evolves the radiation subsystem within the DualFLD module.  
// This is performed in a robust manner; if the current radiation subsystem 
// cannot be solved to the desired tolerance (meaning the time step is 
// likely much too large), it will return a flag to the calling routine to 
// have the step recomputed with a smaller time step size.
int DualFLD::RadStep(HierarchyEntry *ThisGrid, int eta_set, int UV_Xr) {

  // initialize return value
  int recompute_step = 0;

#ifdef USE_PAT
  int rhd_tag;
  char *rhd_label;
  rhd_label = new char[MAX_LINE_LENGTH];
  strcpy(rhd_label, "Hypre_Detail");
  rhd_tag = 1;
#endif /* USE_PAT */

  // update internal Enzo units for current times
  float TempUnits, RadUnits, MassUnits;
  DenUnits0 = LenUnits0 = TempUnits = TimeUnits = VelUnits = MassUnits = 1.0;
  if (GetUnits(&DenUnits0, &LenUnits0, &TempUnits, &TimeUnits, 
	       &VelUnits, &MassUnits, told) != SUCCESS) {
    fprintf(stderr,"DualFLD_RadStep: Error in GetUnits.\n");
    return FAIL; }
  if (RadiationGetUnits(&RadUnits, told) != SUCCESS) {
    fprintf(stderr,"DualFLD_RadStep: Error in RadiationGetUnits.\n");
    return FAIL; }
  float mp = 1.67262171e-24;   // Mass of a proton [g]
  UVUnits0 = RadUnits*UVScale;
  XrUnits0 = RadUnits*XrScale;
  NiUnits0 = (Nchem == 0) ? 1.0 : DenUnits0/mp;
  if (ComovingCoordinates) 
    if (CosmologyComputeExpansionFactor(told, &a0, &adot0) != SUCCESS) {
      fprintf(stderr,"DualFLD_RadStep: Error in CosmologyComputeExpansionFactor.\n");
      return FAIL; }
  adot0 /= TimeUnits;  // rescale to physical units

  DenUnits = LenUnits = TempUnits = TimeUnits = VelUnits = MassUnits = 1.0;
  if (GetUnits(&DenUnits, &LenUnits, &TempUnits, &TimeUnits, 
	       &VelUnits, &MassUnits, tnew) != SUCCESS) {
    fprintf(stderr,"DualFLD_RadStep: Error in GetUnits.\n");
    return FAIL; }
  if (RadiationGetUnits(&RadUnits, tnew) != SUCCESS) {
    fprintf(stderr,"DualFLD_RadStep: Error in RadiationGetUnits.\n");
    return FAIL; }
  UVUnits = RadUnits*UVScale;
  XrUnits = RadUnits*XrScale;
  NiUnits = (Nchem == 0) ? 1.0 : DenUnits/mp;
  if (ComovingCoordinates) 
    if (CosmologyComputeExpansionFactor(tnew, &a, &adot) != SUCCESS) {
      fprintf(stderr,"DualFLD_RadStep: Error in CosmologyComputeExpansionFactor.\n");
      return FAIL; }
  adot /= TimeUnits;  // rescale to physical units
    
  // rescale dt, told, tnew, adot to physical values for use within solver
  dt   *= TimeUnits;
  told *= TimeUnits;
  tnew *= TimeUnits;

  //   compute emissivity at this internal time step (if provided internally)
  if (eta_set == 0) {
    if (this->RadiationSource() != SUCCESS) {
      fprintf(stderr,"DualFLD_RadStep: Error in RadiationSource routine\n");
      return FAIL; }
#ifdef SOLUTION_DIAGNOSTICS
    float Xr_norm=0.0, UV_norm=0.0;
    if (UV_Xr == 0) {
      Xr_norm = extsrc->rmsnorm_component(0);
      if (debug) 
	printf("  Xray emissivity norm = %g\n", Xr_norm);
    } else {
      UV_norm = extsrc->rmsnorm_component(1);
      if (debug) 
	printf("  UV emissivity norm = %g\n",   UV_norm);
    }
#endif
  }
    
  // Multigrid solver: for periodic dims, only coarsen until grid not divisible by 2
  Eint32 max_levels, level=-1;
  int Ndir;
  if (XrBdryType[0][0] == 0) {
    level = 0;
    Ndir = GlobDims[0];
    while ( Ndir%2 == 0 ) {
      level++;
      Ndir /= 2;
    }
  }
  max_levels = level;
  if (rank > 1) {
    if (XrBdryType[1][0] == 0) {
      level = 0;
      Ndir = GlobDims[0];
      while ( Ndir%2 == 0 ) {
	level++;
	Ndir /= 2;
      }
    }
    max_levels = min(level,max_levels);
  }
  if (rank > 2) {
    if (XrBdryType[2][0] == 0) {
      level = 0;
      Ndir = GlobDims[0];
      while ( Ndir%2 == 0 ) {
	level++;
	Ndir /= 2;
      }
    }
    max_levels = min(level,max_levels);
  }
  //  printf("max_levels = %i\n", max_levels);


  //   enforce boundary conditions on current time step vector
  if (this->EnforceBoundary(U0) != SUCCESS) {
    fprintf(stderr,"DualFLD_RadStep: EnforceBoundary failure!!\n");
    return FAIL; }
    

  //   set initial guess as old solution
  sol->copy_component(U0, UV_Xr);
    
  // access relevant chemical species (some may be NULL)
  float *HI   = ThisGrid->GridData->AccessHIDensity();
  float *HeI  = ThisGrid->GridData->AccessHeIDensity();
  float *HeII = ThisGrid->GridData->AccessHeIIDensity();
  //    check that we accessed the required species for Nchem
  if (Nchem > 0) 
    if (HI == NULL) {
      fprintf(stderr,"DualFLD Evolve: cannot access HI density!\n");
      return FAIL; }
  if (Nchem == 3) {
    if (HeI == NULL) {
      fprintf(stderr,"DualFLD Evolve: cannot access HeI density!\n");
      return FAIL;  }
    if (HeII == NULL) {
      fprintf(stderr,"DualFLD Evolve: cannot access HeII density!\n");
      return FAIL; }
  }

  //   compute updated opacities
  if (this->ComputeOpacity(HI, HeI, HeII) != SUCCESS) {
    fprintf(stderr,"DualFLD_RadStep: Error in Opacity routine\n");
    return FAIL; }
    

  // HYPRE stuff
 
  // set up and solve radiation equation
  float *Enew = sol->GetData(UV_Xr);    // updated radiation energy array
  float rhsnorm;                        // used for setting HYPRE solver tolerance
  if (this->SetupSystem(&rhsnorm, UV_Xr) != SUCCESS) {
    fprintf(stderr,"DualFLD_RadStep: Error in SetupSystem routine\n");
    return FAIL; }

  // set solver parameters based on Xray vs UV physics
  float sol_tolerance;
  Eint32 sol_MGmaxit, sol_PCGmaxit, sol_rlxtype, sol_npre, sol_npost;
  if (UV_Xr == 0) {   // Xray
    sol_tolerance = sol_tolerance_Xr;
    sol_MGmaxit   = sol_MGmaxit_Xr;
    sol_PCGmaxit  = sol_PCGmaxit_Xr;
    sol_rlxtype   = sol_rlxtype_Xr;
    sol_npre      = sol_npre_Xr;
    sol_npost     = sol_npost_Xr;
  } else {            // UV
    sol_tolerance = sol_tolerance_UV;
    sol_MGmaxit   = sol_MGmaxit_UV;
    sol_PCGmaxit  = sol_PCGmaxit_UV;
    sol_rlxtype   = sol_rlxtype_UV;
    sol_npre      = sol_npre_UV;
    sol_npost     = sol_npost_UV;
  }


  // only proceed with solve if rhsnorm > eps, otherwise solution increment 
  // is zero and the "solution" to the step is the initial guess
  if (rhsnorm > 1.e-12) {

    // assemble matrix
    Eint32 entries[7] = {0, 1, 2, 3, 4, 5, 6};   // matrix stencil entries
    Eint32 ilower[3]  = {SolvIndices[0][0],SolvIndices[1][0],SolvIndices[2][0]};
    Eint32 iupper[3]  = {SolvIndices[0][1],SolvIndices[1][1],SolvIndices[2][1]};
    
    Eint32 pat_tag = 1;
    Eint32 pat_res = -1;
    
#ifdef USE_PAT
    pat_res = PAT_region_begin(pat_tag, "Hypre_Detail");
    fprintf(stderr, "Initialise_Hypre_Detail %d\n", pat_res);
    // PAT_record(PAT_STATE_ON);
#endif /* USE_PAT */
    
    HYPRE_StructMatrixSetBoxValues(P, ilower, iupper, stSize, entries, matentries); 
    HYPRE_StructMatrixAssemble(P);
    
    // insert rhs into HYPRE vector b
    HYPRE_StructVectorSetBoxValues(rhsvec, ilower, iupper, rhsentries);
    
    // set linear solver tolerance (rescale to relative residual and not actual)
    int use_abs_resid = (rhsnorm < 1.e-8);
    Eflt64 delta = (use_abs_resid) ? sol_tolerance : sol_tolerance/rhsnorm;
    
    // insert sol initial guess into HYPRE vector x 
    int xBuff, yBuff, zBuff, Zbl, Ybl, ix, iy, iz;  // mesh indexing shortcuts
    ilower[0] = SolvIndices[0][0];
    iupper[0] = SolvIndices[0][1];
    xBuff = GhDims[0][0]-SolvOff[0];
    yBuff = (GhDims[1][0]-SolvOff[1])-SolvIndices[1][0];
    zBuff = (GhDims[2][0]-SolvOff[2])-SolvIndices[2][0];
    for (iz=SolvIndices[2][0]; iz<=SolvIndices[2][1]; iz++) {
      Zbl = (iz+zBuff)*ArrDims[0]*ArrDims[1];  ilower[2] = iz;  iupper[2] = iz;
      for (iy=SolvIndices[1][0]; iy<=SolvIndices[1][1]; iy++) {
	Ybl = (iy+yBuff)*ArrDims[0];  ilower[1] = iy;  iupper[1] = iy;
	for (ix=0; ix<=SolvIndices[0][1]-SolvIndices[0][0]; ix++) 
	  HYPREbuff[ix] = 0.0;
	HYPRE_StructVectorSetBoxValues(solvec, ilower, iupper, HYPREbuff);
      }
    }
    
    // assemble vectors
    HYPRE_StructVectorAssemble(solvec);
    HYPRE_StructVectorAssemble(rhsvec);

    // set up the preconditioner
    HYPRE_StructSolver precond;
    HYPRE_StructPFMGCreate(MPI_COMM_WORLD, &precond);
    if (max_levels > -1) 
      HYPRE_StructPFMGSetMaxLevels(precond,  max_levels);
    HYPRE_StructPFMGSetMaxIter(precond,      sol_MGmaxit);
    HYPRE_StructPFMGSetRelaxType(precond,    sol_rlxtype);
    HYPRE_StructPFMGSetNumPreRelax(precond,  sol_npre);
    HYPRE_StructPFMGSetNumPostRelax(precond, sol_npost);

    // set up the solver
    HYPRE_StructSolver solver;
    HYPRE_StructPCGCreate(MPI_COMM_WORLD,  &solver);
    HYPRE_StructPCGSetMaxIter(solver, sol_PCGmaxit);
//    HYPRE_StructPCGSetRelChange(solver, 1);
    HYPRE_StructPCGSetPrintLevel(solver, sol_printl);
    HYPRE_StructPCGSetLogging(solver, sol_log);

    //    set dimension-specific solver options
    if (rank > 1) {
      HYPRE_StructPCGSetPrecond(solver, 
    				(HYPRE_PtrToStructSolverFcn) HYPRE_StructPFMGSolve,  
    				(HYPRE_PtrToStructSolverFcn) HYPRE_StructPFMGSetup, 
    				precond);
    } else {    // ignore preconditioner for 1D tests (bug); increase CG its
      HYPRE_StructPCGSetMaxIter(solver, sol_PCGmaxit*500);
    }
    
    //    set tolerance-specific solver options
    if (use_abs_resid) {
      HYPRE_StructPCGSetTol(solver, 0.0);
      HYPRE_StructPCGSetAbsoluteTol(solver, delta);
    } else {
      HYPRE_StructPCGSetTol(solver, delta);
    }

    // set up and solve the linear system
    HYPRE_StructPCGSetup(solver, P, rhsvec, solvec);
    HYPRE_StructPCGSolve(solver, P, rhsvec, solvec);

    // extract solver & preconditioner statistics
    Eflt64 finalresid=1.0;  // HYPRE solver statistics
    Eint32 Sits=0, Pits=0;  // HYPRE solver statistics
    if (use_abs_resid) 
      finalresid = 0.0;
    else
      HYPRE_StructPCGGetFinalRelativeResidualNorm(solver, &finalresid);
    HYPRE_StructPCGGetNumIterations(solver,   &Sits);
    HYPRE_StructPFMGGetNumIterations(precond, &Pits);
    totIters += Sits;
    if (debug) printf("   lin resid = %.1e (tol = %.1e, |rhs| = %.1e), its = (%i,%i)\n",
		      finalresid*rhsnorm, delta, rhsnorm, Sits, Pits);
    // check if residual is NaN
    if (finalresid != finalresid) {

      // can the step size shrink?
      if (dt > mindt*TimeUnits*1.00000001) {
	// allow remainder of function to complete (to reset units, etc.), 
	// but have calling routine update dt and compute step again.
	recompute_step = 1;
      } else {
	fprintf(stderr,"DualFLD_RadStep: could not solve system at minimum step size!\n");
	
	// output linear system to disk
	if (debug)  fprintf(stderr,"Writing out matrix to file P.mat\n");
	HYPRE_StructMatrixPrint("P.mat",P,0);
	if (debug)  fprintf(stderr,"Writing out rhs to file b.vec\n");
	HYPRE_StructVectorPrint("b.vec",rhsvec,0);
	if (debug)  fprintf(stderr,"Writing out current solution to file x.vec\n");
	HYPRE_StructVectorPrint("x.vec",solvec,0);
	
	fprintf(stderr, "Error in DualFLD_RadStep\n");
	return FAIL;
      }
    }

    // extract values from solution vector, adding them to solution
    for (iz=SolvIndices[2][0]; iz<=SolvIndices[2][1]; iz++) {
      Zbl = (iz+zBuff)*ArrDims[0]*ArrDims[1];
      ilower[2] = iz;  iupper[2] = iz;
      for (iy=SolvIndices[1][0]; iy<=SolvIndices[1][1]; iy++) {
	Ybl = (iy+yBuff)*ArrDims[0];
	ilower[1] = iy;  iupper [1] = iy;
	HYPRE_StructVectorGetBoxValues(solvec, ilower, iupper, HYPREbuff);
	for (ix=0; ix<=SolvIndices[0][1]-SolvIndices[0][0]; ix++) {
	  Enew[Zbl+Ybl+xBuff+ix] += HYPREbuff[ix];
	}
      }
    }
  
    // destroy HYPRE solver & preconditioner structures
    HYPRE_StructPCGDestroy(solver);
    HYPRE_StructPFMGDestroy(precond);

    
#ifdef USE_PAT
    // PAT_record(PAT_STATE_OFF);
    pat_res = PAT_region_end(pat_tag);
    fprintf(stderr, "Finalise_Hypre_Detail %d\n", pat_res);
#endif /* USE_PAT */
    
  } else {
    if (debug) printf("   solve not needed since guess sufficient\n");
  }

  // enforce a solution floor on radiation
  float epsilon=1.0;      // radiation floor
  while (epsilon*0.25 > 0.0)  epsilon*=0.5;
  for (int i=0; i<ArrDims[0]*ArrDims[1]*ArrDims[2]; i++)  
    Enew[i] = max(Enew[i], epsilon);

  // rescale dt, told, tnew, adot back to normalized values
  dt   /= TimeUnits;
  told /= TimeUnits;
  tnew /= TimeUnits;

  return recompute_step;
}

#endif
