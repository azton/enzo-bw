/*****************************************************************************
 * Copyright 2013 Daniel R. Reynolds                                         *
 *                                                                           *
 * This software is released under the terms of the "Enzo Public License"    *
 * in the accompanying LICENSE file.                                         *
 *****************************************************************************/
/***********************************************************************
   Two-Group (X-ray + UV) Flux-Limited Diffusion Implicit Problem Class
   Initialization routine

   author: Daniel R. Reynolds
   date:   January 2013

   PURPOSE: Allocates all necessary internal memory for problem 
            definition and associated linear solver.  This begins the
            interface between Enzo and the FLD solver module, so any
            and all grid/index transformations must be performed and 
            stored here.
************************************************************************/
#ifdef RAD_HYDRO
#include "DualFLD.h"
// #ifdef _OPENMP
// #include "omp.h"
// #endif

// function prototypes
int InitializeRateData(FLOAT Time);
int DualCosmoIonizationInitialize(FILE *fptr, FILE *Outfptr,
				  HierarchyEntry &TopGrid,
				  TopGridData &MetaData, int local);
int DualRHIonizationTestInitialize(FILE *fptr, FILE *Outfptr,
				   HierarchyEntry &TopGrid,
				   TopGridData &MetaData, int local);
int DualRadStreamTestInitialize(FILE *fptr, FILE *Outfptr,
				HierarchyEntry &TopGrid,
				TopGridData &MetaData, int local);
int DualRadConstTestInitialize(FILE *fptr, FILE *Outfptr,
			       HierarchyEntry &TopGrid,
			       TopGridData &MetaData, int local);


int DualFLD::Initialize(HierarchyEntry &TopGrid, TopGridData &MetaData) {

  if (debug)  printf("Entering DualFLD::Initialize routine\n");

  // find the grid corresponding to this process from the Hierarcy
  HierarchyEntry *ThisGrid = &TopGrid;
  int i, dim, face, foundgrid=0;
  for (i=0; i<=MAX_NUMBER_OF_SUBGRIDS; i++) {
    if (MyProcessorNumber != ThisGrid->GridData->ReturnProcessorNumber()) 
      ThisGrid = ThisGrid->NextGridThisLevel;
    else {foundgrid=1; break;}
  }
  if (foundgrid == 0) {
    printf("FLD Initialize ERROR: p%"ISYM" could not locate his grid\n",
	   MyProcessorNumber);
    return FAIL;
  }

// #ifdef _OPENMP
//   // output number of OpenMP threads that can be used in this run
//   int nthreads = omp_get_max_threads();
//   if (debug)
//     printf("DualFLD Initialize: tasks have %"ISYM" available OpenMP threads\n", 
// 	   nthreads);
// #endif  

  // set rank of self-gravity problem to 3
  rank = MetaData.TopGridRank;

  // get processor layout from Grid
  for (dim=0; dim<rank; dim++) 
    layout[dim] = ThisGrid->GridData->GetProcessorLayout(dim);
  
  // get processor location in MPI grid
  for (dim=0; dim<rank; dim++) 
    location[dim] = ThisGrid->GridData->GetProcessorLocation(dim);

  // get neighbor information from grid
  for (dim=0; dim<rank; dim++) 
    for (face=0; face<2; face++) 
      NBors[dim][face] = ThisGrid->GridData->GetProcessorNeighbors(dim,face);

  // set default module parameters
  int XrOnly   = 0;    // UV + Xray propagation
  int staticXr = 0;    // dynamic Xray radiation
  int staticUV = 0;    // dynamic UV radiation
  Nchem       = 3;     // hydrogen + helium
  UVSpectrum  = -1;    // monochromatic spectrum
  UVFrequency = 13.6;  // monochromatic spectrum frequency (eV)
  XrSpectrum  = -1;    // monochromatic spectrum
  XrFrequency = 500.0; // monochromatic spectrum frequency (500 keV)
  HFrac       = 0.76;  // all Hydrogen
  theta       = 1.0;   // backwards euler implicit time discret.
  dtnorm      = 2.0;   // use 2-norm for time step estimation
  dtgrowth    = 1.1;   // 10% allowed growth in dt per step
  UVScale     = 1.0;   // no radiation equation scaling
  XrScale     = 1.0;   // no radiation equation scaling
  int autoscale = 1;   // enable automatic variable scaling
  for (dim=0; dim<rank; dim++)     // set default radiation boundaries to 
    for (face=0; face<2; face++) { //   periodic in each direction
      XrBdryType[dim][face] = 0;
      UVBdryType[dim][face] = 0;
    }

  // set default solver parameters
  sol_tolerance_Xr   = 1.0e-8;    // HYPRE solver tolerance
  sol_tolerance_UV   = 1.0e-8;
  sol_MGmaxit_Xr     = 5;         // HYPRE max multigrid iters
  sol_MGmaxit_UV     = 3;
  sol_PCGmaxit_Xr    = 3;         // HYPRE max PCG iters
  sol_PCGmaxit_UV    = 2;
  sol_rlxtype_Xr     = 2;         // HYPRE relaxation type
  sol_rlxtype_UV     = 1;
  sol_npre_Xr        = 3;         // HYPRE num pre-smoothing steps
  sol_npre_UV        = 1;
  sol_npost_Xr       = 3;         // HYPRE num post-smoothing steps
  sol_npost_UV       = 1;
  sol_printl         = 1;         // HYPRE print level
  sol_log            = 1;         // HYPRE logging level

  // set default ionization parameters
  NGammaDotUV        = 0.0;       // UV source strength
  NGammaDotXr        = 0.0;       // X-ray source strength
  EtaRadius          = 0.0;       // single cell
  EtaCenter[0]       = 0.0;       // x-location
  EtaCenter[1]       = 0.0;       // y-location
  EtaCenter[2]       = 0.0;       // z-location

  // set default chemistry constants
  hnu0_HI   = 13.6;      // ionization energy of HI   [eV]
  hnu0_HeI  = 24.6;      // ionization energy of HeI  [eV]
  hnu0_HeII = 54.4;      // ionization energy of HeII [eV]


  ////////////////////////////////
  // if input file present, over-write defaults with module inputs
  FILE *fptr;
  char line[MAX_LINE_LENGTH];
  int ret;
  char *dummy = new char[MAX_LINE_LENGTH];
  dummy[0] = 0;

  // check whether input file is non-null
  if (MetaData.RadHydroParameterFname != NULL) {
    if ((fptr = fopen(MetaData.RadHydroParameterFname, "r")) == NULL)
      fprintf(stderr,"Error opening RadHydro parameter file %s, using defaults\n",
	      MetaData.RadHydroParameterFname);
    else {

      // read until out of lines
      while (fgets(line, MAX_LINE_LENGTH, fptr) != NULL) {
	ret = 0;
	ret += sscanf(line, "DualFLDXrayOnly = %"ISYM, &XrOnly);
	ret += sscanf(line, "DualFLDXrayStatic = %"ISYM, &staticXr);
	ret += sscanf(line, "DualFLDUVStatic = %"ISYM, &staticUV);
	ret += sscanf(line, "DualFLDUVSpectrum = %"ISYM, &UVSpectrum);
	ret += sscanf(line, "DualFLDUVFrequency = %"FSYM, &UVFrequency);
	ret += sscanf(line, "DualFLDXraySpectrum = %"ISYM, &XrSpectrum);
	ret += sscanf(line, "DualFLDXrayFrequency = %"FSYM, &XrFrequency);
	ret += sscanf(line, "DualFLDChemistry = %"ISYM, &Nchem);
	ret += sscanf(line, "DualFLDHFraction = %"FSYM, &HFrac);
	ret += sscanf(line, "DualFLDMaxDt = %"FSYM, &maxdt);
	ret += sscanf(line, "DualFLDMinDt = %"FSYM, &mindt);
	ret += sscanf(line, "DualFLDInitDt = %"FSYM, &initdt);
	ret += sscanf(line, "DualFLDDtNorm = %"FSYM, &dtnorm);
	ret += sscanf(line, "DualFLDDtGrowth = %"FSYM, &dtgrowth);
	ret += sscanf(line, "DualFLDDtUVFac = %"FSYM, &dtfac[0]);
	ret += sscanf(line, "DualFLDDtXrayFac = %"FSYM, &dtfac[1]);
	ret += sscanf(line, "DualFLDUVScaling = %"FSYM, &UVScale);
	ret += sscanf(line, "DualFLDXrayScaling = %"FSYM, &XrScale);
	ret += sscanf(line, "DualFLDAutoScaling = %"ISYM, &autoscale);
	ret += sscanf(line, "DualFLDTheta = %"FSYM, &theta);
	ret += sscanf(line, "DualFLDUVBoundaryX0Faces = %"ISYM" %"ISYM, 
		      UVBdryType[0], UVBdryType[0]+1);
	if (rank > 1) {
	  ret += sscanf(line, "DualFLDUVBoundaryX1Faces = %"ISYM" %"ISYM,
			UVBdryType[1], UVBdryType[1]+1);
	  if (rank > 2) {
	    ret += sscanf(line, "DualFLDUVBoundaryX2Faces = %"ISYM" %"ISYM,
			  UVBdryType[2], UVBdryType[2]+1);
	  }
	}
	ret += sscanf(line, "DualFLDXrBoundaryX0Faces = %"ISYM" %"ISYM, 
		      XrBdryType[0], XrBdryType[0]+1);
	if (rank > 1) {
	  ret += sscanf(line, "DualFLDXrBoundaryX1Faces = %"ISYM" %"ISYM,
			XrBdryType[1], XrBdryType[1]+1);
	  if (rank > 2) {
	    ret += sscanf(line, "DualFLDXrBoundaryX2Faces = %"ISYM" %"ISYM,
			  XrBdryType[2], XrBdryType[2]+1);
	  }
	}
	ret += sscanf(line, "DualFLDSolToleranceXray = %"FSYM, &sol_tolerance_Xr);
	ret += sscanf(line, "DualFLDSolToleranceUV = %"FSYM, &sol_tolerance_UV);
	ret += sscanf(line, "DualFLDMaxMGItersXray = %i", &sol_MGmaxit_Xr);
	ret += sscanf(line, "DualFLDMaxMGItersUV = %i", &sol_MGmaxit_UV);
	ret += sscanf(line, "DualFLDMaxPCGItersXray = %i", &sol_PCGmaxit_Xr);
	ret += sscanf(line, "DualFLDMaxPCGItersUV = %i", &sol_PCGmaxit_UV);
	ret += sscanf(line, "DualFLDMGRelaxTypeXray = %i", &sol_rlxtype_Xr);
	ret += sscanf(line, "DualFLDMGRelaxTypeUV = %i", &sol_rlxtype_UV);
	ret += sscanf(line, "DualFLDMGPreRelaxXray = %i", &sol_npre_Xr);
	ret += sscanf(line, "DualFLDMGPreRelaxUV = %i", &sol_npre_UV);
	ret += sscanf(line, "DualFLDMGPostRelaxXray = %i", &sol_npost_Xr);
	ret += sscanf(line, "DualFLDMGPostRelaxUV = %i", &sol_npost_UV);
	ret += sscanf(line, "DualFLDNGammaDotXray = %"FSYM, &NGammaDotXr);
	ret += sscanf(line, "DualFLDNGammaDotUV = %"FSYM, &NGammaDotUV);
	ret += sscanf(line, "DualFLDEtaRadius = %"FSYM, &EtaRadius);
	ret += sscanf(line, "DualFLDEtaCenter = %"FSYM" %"FSYM" %"FSYM, 
		      &(EtaCenter[0]), &(EtaCenter[1]), &(EtaCenter[2]));
	
      }  // end loop over file lines
    }  // end successful file open
  }  // end if file name exists
 
  // clean up
  delete[] dummy;
  rewind(fptr);
  fclose(fptr);

  // set XrayOnly flag based on input value XrOnly
  XrayOnly = (XrOnly == 1);

  // set static radiation fields based on input values staticXr and staticUV
  XrStatic = (staticXr == 1);
  UVStatic = (staticUV == 1);


  ////////////////////////////////

  // check that these give appropriate values, otherwise set dim to periodic
  for (dim=0; dim<rank; dim++) 
    for (face=0; face<2; face++) {
      /// ADD NEW BOUNDARY CONDITION TYPES HERE!
      if ((XrBdryType[dim][face] < 0) || (XrBdryType[dim][face] > 2)) {
	fprintf(stderr,"DualFLD_Initialize Warning: re-setting Xray BC to periodic, dim %"ISYM", face %"ISYM"\n",dim,face);
	XrBdryType[dim][face] = 0;
      }
      if (!XrayOnly)
	if ((UVBdryType[dim][face] < 0) || (UVBdryType[dim][face] > 2)) {
	  fprintf(stderr,"DualFLD_Initialize Warning: re-setting UV BC to periodic, dim %"ISYM", face %"ISYM"\n",dim,face);
	  UVBdryType[dim][face] = 0;
	}
    }

  // check that periodic faces match
  for (dim=0; dim<rank; dim++) {
    if ((XrBdryType[dim][0]*XrBdryType[dim][1] == 0) && 
	(XrBdryType[dim][0]+XrBdryType[dim][1] != 0)) {
      fprintf(stderr,"DualFLD_Initialize Warning: non-matching periodic Xray BCs, dim %"ISYM"\n",dim);
      XrBdryType[dim][0] = 0;
      XrBdryType[dim][1] = 0;
    }
    if (!XrayOnly)
      if ((UVBdryType[dim][0]*UVBdryType[dim][1] == 0) && 
	  (UVBdryType[dim][0]+UVBdryType[dim][1] != 0)) {
	fprintf(stderr,"DualFLD_Initialize Warning: non-matching periodic UV BCs, dim %"ISYM"\n",dim);
	UVBdryType[dim][0] = 0;
	UVBdryType[dim][1] = 0;
      }
  }

  // check that Xray and UV periodicity match
  if (!XrayOnly)
    for (dim=0; dim<rank; dim++) {
      if (((XrBdryType[dim][0]*XrBdryType[dim][1] == 0) && 
	   (UVBdryType[dim][0]*UVBdryType[dim][1] != 0)) ||
	  ((UVBdryType[dim][0]*UVBdryType[dim][1] == 0) && 
	   (XrBdryType[dim][0]*XrBdryType[dim][1] != 0))) {
	fprintf(stderr,"DualFLD_Initialize Warning: non-matching periodicity of Xray and UV fields, dim %"ISYM", Xr = (%"ISYM",%"ISYM"), UV = (%"ISYM",%"ISYM")\n", dim, XrBdryType[dim][0], XrBdryType[dim][1], UVBdryType[dim][0], UVBdryType[dim][1]);
	XrBdryType[dim][0] = 0;
	XrBdryType[dim][1] = 0;
	UVBdryType[dim][0] = 0;
	UVBdryType[dim][1] = 0;
      }
    }

  // ensure that new BdryVals array pointers are set to NULL
  for (dim=0; dim<rank; dim++) 
    for (face=0; face<2; face++) {
      if (XrBdryVals[dim][face] != NULL) 
	delete[] XrBdryVals[dim][face];
      XrBdryVals[dim][face] = NULL;
      if (UVBdryVals[dim][face] != NULL) 
	delete[] UVBdryVals[dim][face];
      UVBdryVals[dim][face] = NULL;
    }

  // set up subdomain information
  //   EdgeVals gives the location of the left/right edge of the
  //      domain (no bdry) -- start with Enzo grid size
  for (dim=0; dim<rank; dim++) {
    EdgeVals[dim][0] = ThisGrid->GridData->GetGridLeftEdge(dim);
    EdgeVals[dim][1] = ThisGrid->GridData->GetGridRightEdge(dim);
  }

  //   LocDims holds the dimensions of the local domain, 
  //   active cells only (no ghost or boundary cells)
  for (dim=0; dim<rank; dim++)
    LocDims[dim] = ThisGrid->GridData->GetGridEndIndex(dim)
      - ThisGrid->GridData->GetGridStartIndex(dim) + 1;

  // Nchem gives the number of chemical species
  if ((Nchem < 1) || (Nchem > 3)) {
    fprintf(stderr,"DualFLD Initialize: illegal Nchem = %"ISYM"\n",Nchem);
    fprintf(stderr,"   re-setting Nchem to 3\n");
    Nchem = 3;  // default is hydrogen + helium
  }

  // HFrac must be between 0 and 1
  if ((HFrac < 0.0) || (HFrac > 1.0)) {
    fprintf(stderr,"DualFLD Initialize: illegal DualFLDHFraction = %g\n",
	    HFrac);
    fprintf(stderr,"   re-setting to 1.0\n");
    HFrac = 1.0;  // default is all Hydrogen
  }

  // monochromatic radiation frequencies must be >0
  if ((UVSpectrum < 0) && (UVFrequency <= 0.0)) {
    fprintf(stderr,"DualFLD Initialize: illegal UVFrequency = %g\n",UVFrequency);
    fprintf(stderr,"   re-setting UVFrequency 13.6\n");
    UVFrequency = 13.6;   // default is the hydrogen ionization threshold
  }
  if ((XrSpectrum < 0) && (XrFrequency <= 0.0)) {
    fprintf(stderr,"DualFLD Initialize: illegal XrFrequency = %g\n",XrFrequency);
    fprintf(stderr,"   re-setting XrFrequency 500.0\n");
    XrFrequency = 500.0;  // default is 0.5 keV
  }

  // maxdt gives the maximum radiation time step size
  if (maxdt <= 0.0) {
    fprintf(stderr,"DualFLD Initialize: illegal DualFLDMaxDt = %g\n",maxdt);
    fprintf(stderr,"   re-setting to %g\n",huge_number);
    maxdt = huge_number;  // default is no limit
  }

  // mindt gives the minimum radiation time step size
  if (mindt < 0.0) {
    fprintf(stderr,"DualFLD Initialize: illegal DualFLDMinDt = %g\n",mindt);
    fprintf(stderr,"   re-setting to %g\n",0.0);
    mindt = 0.0;  // default is 0.0
  }

  // initdt gives the initial time step size
  if (initdt <= 0.0) {
    fprintf(stderr,"DualFLD Initialize: illegal DualFLDInitDt = %g\n",initdt);
    fprintf(stderr,"   re-setting to %g\n",huge_number);
    initdt = huge_number;  // default is no limit
  }

  // a, adot give cosmological expansion & rate
  a = 1.0;  a0 = 1.0;  adot = 0.0;  adot0 = 0.0;

  // *Scale give variable scalings for implicit solver
  if (UVScale <= 0.0) {
    fprintf(stderr,"Initialize: illegal DualFLDUVScaling = %g\n",UVScale);
    fprintf(stderr,"   re-setting to 1.0\n");
    UVScale = 1.0;  // default is no scaling
  }
  if (XrScale <= 0.0) {
    fprintf(stderr,"Initialize: illegal DualFLDXrayScaling = %g\n",XrScale);
    fprintf(stderr,"   re-setting to 1.0\n");
    XrScale = 1.0;  // default is no scaling
  }
  autoScale = (autoscale != 0);  // set bool based on integer input
  if (debug)
    printf("DualFLD::Initialize p%"ISYM": UVScale = %g, XrScale = %g, autoScale = %i\n",
	   MyProcessorNumber,UVScale,XrScale,autoscale);

  // dtfac gives the desired percent change in values per step
  if (dtfac[0] <= 0.0) {
    fprintf(stderr,"DualFLD Initialize: illegal DualFLDDtUVFac = %g\n",dtfac[0]);
    fprintf(stderr,"   re-setting to %g\n",huge_number);
    dtfac[0] = huge_number;  // default is no limit
  }
  if (dtfac[1] <= 0.0) {
    fprintf(stderr,"DualFLD Initialize: illegal DualFLDDtXrayFac = %g\n",dtfac[1]);
    fprintf(stderr,"   re-setting to %g\n",huge_number);
    dtfac[1] = huge_number;  // default is no limit
  }

  // dtnorm gives the norm for calculating allowed relative change per step
  if (dtnorm < 0.0) {
    fprintf(stderr,"DualFLD Initialize: illegal DualFLDDtNorm = %g\n",dtnorm);
    fprintf(stderr,"   re-setting to 2.0 (2-norm)\n");
    dtnorm = 2.0;  // default is 2-norm
  }

  // dtgrowth gives the maximum growth factor in dt per step
  if (dtgrowth < 1.0 || dtgrowth > 10.0) {
    fprintf(stderr,"DualFLD Initialize: illegal DualFLDDtGrowth = %g\n",dtgrowth);
    fprintf(stderr,"   re-setting to 1.1\n");
    dtgrowth = 1.1;
  }

  // theta gives the implicit time-stepping method (1->BE, 0.5->CN, 0->FE)
  if ((theta < 0.0) || (theta > 1.0)) {
    fprintf(stderr,"DualFLD Initialize: illegal DualFLDTheta = %g\n",
	    theta);
    fprintf(stderr,"   re-setting theta to 1.0 (Backwards Euler)\n");
    theta = 1.0;  // default is backwards Euler
  }

  // set flags denoting if this processor is on the external boundary
  for (dim=0; dim<rank; dim++) {
    if (layout[dim]==0) {
      OnBdry[dim][0] = OnBdry[dim][1] = true;
    }
    else {
      OnBdry[dim][0] = (location[dim] == 0);
      OnBdry[dim][1] = (location[dim] == layout[dim]-1);
    }
  }
  if (debug){
    printf("DualFLD::Initialize p%"ISYM": rank = %"ISYM", Nchem = %"ISYM", HFrac = %g\n",
	   MyProcessorNumber, rank, Nchem, HFrac);
    printf("DualFLD::Initialize p%"ISYM": layout = (%"ISYM",%"ISYM",%"ISYM")\n",
	   MyProcessorNumber,layout[0],layout[1],layout[2]);
  }

  //   for non-periodic domain, unset neighbor info.
#ifndef USE_MPI
  int MPI_PROC_NULL = -3;
#endif
  for (dim=0; dim<rank; dim++) {
    if ((OnBdry[dim][0]) && (XrBdryType[dim][0] != 0))
      NBors[dim][0] = MPI_PROC_NULL;
    if ((OnBdry[dim][1]) && (XrBdryType[dim][1] != 0))
      NBors[dim][1] = MPI_PROC_NULL;
  }
  if (debug) {
    printf("DualFLD::Initialize p%"ISYM": XrBdryType = (%"ISYM":%"ISYM",%"ISYM":%"ISYM",%"ISYM":%"ISYM")\n",
	   MyProcessorNumber, XrBdryType[0][0], XrBdryType[0][1], XrBdryType[1][0], 
	   XrBdryType[1][1], XrBdryType[2][0], XrBdryType[2][1]);
    printf("DualFLD::Initialize p%"ISYM": UVBdryType = (%"ISYM":%"ISYM",%"ISYM":%"ISYM",%"ISYM":%"ISYM")\n",
	   MyProcessorNumber, UVBdryType[0][0], UVBdryType[0][1], UVBdryType[1][0], 
	   UVBdryType[1][1], UVBdryType[2][0], UVBdryType[2][1]);
  }

  // initialize the time values
  tnew = told = MetaData.Time;

  // dt* gives the time step sizes for each piece of physics
  dtrad = initdt;             // use the input value (scaled units)

  // set the global initial dt to match the radiation timestep
  dt = initdt;
  ThisGrid->GridData->SetMaxRadiationDt(dt);
  
  // compute global dimension information
  for (dim=0; dim<rank; dim++)
    GlobDims[dim] = MetaData.TopGridDims[dim];

  // dx gives grid cell size (comoving, normalized units)
  for (dim=0; dim<rank; dim++)
    dx[dim] = (EdgeVals[dim][1]-EdgeVals[dim][0])/LocDims[dim];

  // compute global index information for this subdomain
  float fCellsLeft;
  for (dim=0; dim<rank; dim++) {

    // the global indexing is easy if we're at the left edge
    if (location[dim]==0)  SolvIndices[dim][0]=0;

    // otherwise we compute the number of intervening cells to left edge
    else {

      // get floating point value for number of cells
      fCellsLeft = (EdgeVals[dim][0] - DomainLeftEdge[dim])/dx[dim];

      // round floating point value to closest integer
      SolvIndices[dim][0] =  (long) (fCellsLeft >= 0.0) ?
	(trunc(fCellsLeft+0.5)) : (trunc(fCellsLeft-0.5));
    }

    // add on local size to obtain right edge indices
    SolvIndices[dim][1] = SolvIndices[dim][0] + LocDims[dim] - 1;
  }

  // store local array sizes (active + ghost)
  for (dim=0; dim<rank; dim++)
    ArrDims[dim] = LocDims[dim] + 2*DEFAULT_GHOST_ZONES;

  // set up vector container for previous time step (empty data)
  int xghosts = DEFAULT_GHOST_ZONES, yghosts=0, zghosts=0;
  if (rank > 1) {
    yghosts = DEFAULT_GHOST_ZONES;
    if (rank > 2) {
      zghosts = DEFAULT_GHOST_ZONES;
    }
  }
  int empty=1;
  int ngroups = (XrayOnly) ? 1 : 2;
  U0 = new EnzoVector(LocDims[0], LocDims[1], LocDims[2], 
		      xghosts, xghosts, yghosts, yghosts, zghosts, zghosts, 
		      ngroups, NBors[0][0], NBors[0][1], NBors[1][0], 
		      NBors[1][1], NBors[2][0], NBors[2][1], empty);
  GhDims[0][0] = xghosts;
  GhDims[0][1] = xghosts;
  GhDims[1][0] = yghosts;
  GhDims[1][1] = yghosts;
  GhDims[2][0] = zghosts;
  GhDims[2][1] = zghosts;

  // set up vectors for temporary storage and Jacobian components
  sol     = U0->clone();
  extsrc  = U0->clone();
  opacity = U0->clone();


  // ensure that CoolData object has been set up, and reset Hydrogen fraction
  if (CoolData.ceHI == NULL) 
    if (InitializeRateData(MetaData.Time) == FAIL) {
      fprintf(stderr,"Error in InitializeRateData.\n");
      return FAIL;
    }
  CoolData.HydrogenFractionByMass = HFrac;


  // compute Radiation Energy spectrum integrals
  if (this->ComputeRadiationIntegrals() == FAIL) {
    fprintf(stderr,"DualFLD::Initialize Error in computing radiation spectrum integrals\n");
    return FAIL;
  }

  // initialize HYPRE stuff
  float stime = MPI_Wtime();

  //    initialize the diagnostic information
  totIters = 0;

  //    set up the grid
  //       create the grid object
  HYPRE_StructGridCreate(MPI_COMM_WORLD, rank, &grid);

  //       set my grid extents as if we have one part with multiple boxes.
  //       Have each processor describe it's own global extents
  Eint32 ilower[3] = {SolvIndices[0][0], SolvIndices[1][0], SolvIndices[2][0]};
  Eint32 iupper[3] = {SolvIndices[0][1], SolvIndices[1][1], SolvIndices[2][1]};
  HYPRE_StructGridSetExtents(grid, ilower, iupper);

  //       set grid periodicity
  Eint32 periodicity[3] = {0, 0, 0};
  if (XrBdryType[0][0] == 0)  periodicity[0] = GlobDims[0];
  if (XrBdryType[1][0] == 0)  periodicity[1] = GlobDims[1];
  if (XrBdryType[2][0] == 0)  periodicity[2] = GlobDims[2];
  HYPRE_StructGridSetPeriodic(grid, periodicity);
  
  //       assemble the grid
  HYPRE_StructGridAssemble(grid);

  //   set up the stencil
  if (rank == 1) 
    stSize = 3;
  else if (rank == 2)
    stSize = 5;
  else 
    stSize = 7;
  HYPRE_StructStencilCreate(rank, stSize, &stencil);

  //      set stencil entries
  Eint32 offset[3];
  Eint32 stentry=0;
  //         dependency to x2 left
  if (rank == 3) {
    offset[0] = 0;  offset[1] = 0;  offset[2] = -1;
    HYPRE_StructStencilSetElement(stencil, stentry++, offset);
  }
  //         dependency to x1 left
  if (rank >= 2) {
    offset[0] = 0;  offset[1] = -1;  offset[2] = 0;
    HYPRE_StructStencilSetElement(stencil, stentry++, offset);
  }
  //         dependency to x0 left
  offset[0] = -1;  offset[1] = 0;  offset[2] = 0;
  HYPRE_StructStencilSetElement(stencil, stentry++, offset);
  //         dependency to self
  offset[0] = 0;  offset[1] = 0;  offset[2] = 0;
  HYPRE_StructStencilSetElement(stencil, stentry++, offset);
  //         dependency to x0 right
  offset[0] = 1;  offset[1] = 0;  offset[2] = 0;
  HYPRE_StructStencilSetElement(stencil, stentry++, offset);
  //         dependency to x1 right
  if (rank >= 2) {
    offset[0] = 0;  offset[1] = 1;  offset[2] = 0;
    HYPRE_StructStencilSetElement(stencil, stentry++, offset);
  }
  //         dependency to x2 right
  if (rank == 3) {
    offset[0] = 0;  offset[1] = 0;  offset[2] = 1;
    HYPRE_StructStencilSetElement(stencil, stentry++, offset);
  }

  //   allocate temporary arrays
  int Nx = (SolvIndices[0][1]-SolvIndices[0][0]+1);
  int Ny = (SolvIndices[1][1]-SolvIndices[1][0]+1);
  int Nz = (SolvIndices[2][1]-SolvIndices[2][0]+1);
  matentries = new Eflt64[stSize*Nx*Ny*Nz];
  rhsentries = new Eflt64[Nx*Ny*Nz];
  HYPREbuff  = new Eflt64[Nx];
  HYPRE_StructMatrixCreate(MPI_COMM_WORLD, grid, stencil, &P);
  HYPRE_StructMatrixInitialize(P);
  HYPRE_StructVectorCreate(MPI_COMM_WORLD, grid, &rhsvec);
  HYPRE_StructVectorInitialize(rhsvec);
  HYPRE_StructVectorCreate(MPI_COMM_WORLD, grid, &solvec);
  HYPRE_StructVectorInitialize(solvec);

  float ftime = MPI_Wtime();
  HYPREtime += ftime-stime;


  //   check MG solver parameters
  if (sol_MGmaxit_Xr < 0) {
    fprintf(stderr,"Illegal DualFLDMaxMGItersXray = %i. Setting to 5\n",
	    sol_MGmaxit_Xr);
    sol_MGmaxit_Xr = 5;
  }
  if (sol_MGmaxit_UV < 0) {
    fprintf(stderr,"Illegal DualFLDMaxMGItersUV = %i. Setting to 3\n",
	    sol_MGmaxit_UV);
    sol_MGmaxit_UV = 3;
  }
  if (sol_PCGmaxit_Xr < 0) {
    fprintf(stderr,"Illegal DualFLDMaxPCGItersXray = %i. Setting to 3\n",
	    sol_PCGmaxit_Xr);
    sol_PCGmaxit_Xr = 3;
  }
  if (sol_PCGmaxit_UV < 0) {
    fprintf(stderr,"Illegal DualFLDMaxPCGItersUV = %i. Setting to 2\n",
	    sol_PCGmaxit_UV);
    sol_PCGmaxit_UV = 2;
  }
  if ((sol_rlxtype_Xr<0) || (sol_rlxtype_Xr>3)) {
    fprintf(stderr,"Illegal DualFLDMGRelaxTypeXray = %i. Setting to 2\n",
	    sol_rlxtype_Xr);
    sol_rlxtype_Xr = 2;
  }
  if ((sol_rlxtype_UV<0) || (sol_rlxtype_UV>3)) {
    fprintf(stderr,"Illegal DualFLDMGRelaxTypeUV = %i. Setting to 1\n",
	    sol_rlxtype_UV);
    sol_rlxtype_UV = 1;
  }
  if (sol_npre_Xr < 1) {
    fprintf(stderr,"Illegal DualFLDMGPreRelaxXray = %i. Setting to 3\n",
	    sol_npre_Xr);
    sol_npre_Xr = 3;
  }
  if (sol_npre_UV < 1) {
    fprintf(stderr,"Illegal DualFLDMGPreRelaxUV = %i. Setting to 1\n",
	    sol_npre_UV);
    sol_npre_UV = 1;
  }
  if (sol_npost_Xr < 1) {
    fprintf(stderr,"Illegal DualFLDMGPostRelaxXray = %i. Setting to 3\n",
	    sol_npost_Xr);
    sol_npost_Xr = 3;
  }
  if (sol_npost_UV < 1) {
    fprintf(stderr,"Illegal DualFLDMGPostRelaxUV = %i. Setting to 1\n",
	    sol_npost_UV);
    sol_npost_UV = 1;
  }
  if ((sol_tolerance_Xr < 1.0e-15) || (sol_tolerance_Xr > 1.0)) {
    fprintf(stderr,"Illegal DualFLDSolToleranceXray = %g. Setting to 1e-4\n",
	    sol_tolerance_Xr);
    sol_tolerance_Xr = 1.0e-4;
  }
  if ((sol_tolerance_UV < 1.0e-15) || (sol_tolerance_UV > 1.0)) {
    fprintf(stderr,"Illegal DualFLDSolToleranceUV = %g. Setting to 1e-4\n",
	    sol_tolerance_UV);
    sol_tolerance_UV = 1.0e-4;
  }


  ////////////////////////////////
  // set up the boundary conditions on the radiation field, 
  // depending on the ProblemType
  float ONE = 0.9e-6;
  float ZERO = 0.0;

  // set boundary conditions based on problem type
  // (default to homogeneous Dirichlet)
  switch (ProblemType) {

  // Streaming test problem: set Dirichlet BC to value of 1.0, 
  // or Neumann BC to value of 0.0; leave Periodic BC alone
  case 201:
    // first call local problem initializer (to allocate/setup local data)
    if (DualRadStreamTestInitialize(fptr, fptr, TopGrid, 
				    MetaData, 1) == FAIL) {
      fprintf(stderr,"Error in DualRadStreamTestInitialize.\n");
      return FAIL;
    }
    
    // set boundary conditions at various faces (Dirichlet has value 1, Neumann has 0)
    for (dim=0; dim<MetaData.TopGridRank; dim++) {
      for (face=0; face<2; face++) {
	if (XrBdryType[dim][face] == 1) {
	  printf("Setting Dirichlet condition in Xr at dim %"ISYM", face %"ISYM"\n",
		 dim, face);
	  if (this->SetupBoundary(dim,face,0,1,&ONE) == FAIL) {
	    fprintf(stderr,"Error setting dim %"ISYM", face %"ISYM" Xray radiation BCs\n",
		    dim, face);
	    return FAIL;
	  }
	}
	else if (XrBdryType[dim][face] == 2) {
	  printf("Setting Neumann condition in Xr at dim %"ISYM", face %"ISYM"\n",
		 dim, face);
	  if (this->SetupBoundary(dim,face,0,1,&ZERO) == FAIL) {
	    fprintf(stderr,"Error setting dim %"ISYM", face %"ISYM" Xray radiation BCs\n",
		    dim, face);
	    return FAIL;
	  }
	}
	if (UVBdryType[dim][face] == 1) {
	  printf("Setting Dirichlet condition in UV at dim %"ISYM", face %"ISYM"\n",
		 dim, face);
	  if (this->SetupBoundary(dim,face,1,1,&ONE) == FAIL) {
	    fprintf(stderr,"Error setting dim %"ISYM", face %"ISYM" UV radiation BCs\n",
		    dim, face);
	    return FAIL;
	  }
	}
	else if (UVBdryType[dim][face] == 2) {
	  printf("Setting Neumann condition in UV at dim %"ISYM", face %"ISYM"\n",
		 dim, face);
	  if (this->SetupBoundary(dim,face,1,1,&ZERO) == FAIL) {
	    fprintf(stderr,"Error setting dim %"ISYM", face %"ISYM" UV radiation BCs\n",
		    dim, face);
	    return FAIL;
	  }
	}
      }
    }
    break;
    
    
  // Ionization test initializer
  case 210:
  case 211:
    // first call local problem initializer (to allocate/setup local data)
    if (DualRHIonizationTestInitialize(NULL, NULL, TopGrid, 
				       MetaData, 1) == FAIL) {
      fprintf(stderr,"Error in DualRHIonizationTestInitialize.\n");
      return FAIL;
    }
    
    // set BC on all non-periodic faces to zero-valued (Dirichlet or Neumann)
    for (dim=0; dim<MetaData.TopGridRank; dim++) {
      for (face=0; face<2; face++) {
	if (XrBdryType[dim][face] != 0)
	  if (this->SetupBoundary(dim,face,0,1,&ZERO) == FAIL) {
	    fprintf(stderr,"Error setting dim %"ISYM", face %"ISYM" Xray radiation BCs\n",
		    dim, face);
	    return FAIL;
	  }
	if (UVBdryType[dim][face] != 0)
	  if (this->SetupBoundary(dim,face,1,1,&ZERO) == FAIL) {
	    fprintf(stderr,"Error setting dim %"ISYM", face %"ISYM" UV radiation BCs\n",
		    dim, face);
	    return FAIL;
	  }
      }
    }
    break;


  // Cosmology test initializer (periodic BCs on all faces)
  case 214:
    // first call local problem initializer (to allocate/setup local data)
    if (DualCosmoIonizationInitialize(NULL, NULL, TopGrid, MetaData, 1) == FAIL) {
      fprintf(stderr,"Error in DualCosmoIonizationInitialize.\n");
      return FAIL;
    }
    
    break;


  // Homogeneous test initializer
  case 216:
    // first call local problem initializer (to allocate/setup local data)
    if (DualRadConstTestInitialize(NULL, NULL, TopGrid, MetaData, 1) == FAIL) {
      fprintf(stderr,"Error in DualRadConstTestInitialize.\n");
      return FAIL;
    }
    
    // set BC on all non-periodic faces to zero-valued (Dirichlet or Neumann)
    for (dim=0; dim<MetaData.TopGridRank; dim++) {
      for (face=0; face<2; face++) {
	if (XrBdryType[dim][face] != 0)
	  if (this->SetupBoundary(dim,face,0,1,&ZERO) == FAIL) {
	    fprintf(stderr,"Error setting dim %"ISYM", face %"ISYM" Xray radiation BCs\n",
		    dim, face);
	    return FAIL;
	  }
	if (UVBdryType[dim][face] != 0)
	  if (this->SetupBoundary(dim,face,1,1,&ZERO) == FAIL) {
	    fprintf(stderr,"Error setting dim %"ISYM", face %"ISYM" UV radiation BCs\n",
		    dim, face);
	    return FAIL;
	  }
      }
    }
    break;


    // Insert new problem intializers here...



  default:
    // set BC on all non-periodic faces to zero-valued (Dirichlet or Neumann)
    for (dim=0; dim<MetaData.TopGridRank; dim++) {
      for (face=0; face<2; face++) {
	if (XrBdryType[dim][face] != 0)
	  if (this->SetupBoundary(dim,face,0,1,&ZERO) == FAIL) {
	    fprintf(stderr,"Error setting dim %"ISYM", face %"ISYM" Xray radiation BCs\n",
		    dim, face);
	    return FAIL;
	  }
	if (UVBdryType[dim][face] != 0)
	  if (this->SetupBoundary(dim,face,1,1,&ZERO) == FAIL) {
	    fprintf(stderr,"Error setting dim %"ISYM", face %"ISYM" UV radiation BCs\n",
		    dim, face);
	    return FAIL;
	  }
      }
    }
  }
  ////////////////////////////////

  // initialize rate exchange fields to zero
  if (RadiativeCooling) {
    float *phHI       = ThisGrid->GridData->AccessKPhHI();
    float *phHeI      = ThisGrid->GridData->AccessKPhHeI();
    float *phHeII     = ThisGrid->GridData->AccessKPhHeII();
    float *photogamma = ThisGrid->GridData->AccessPhotoGamma();
    float *dissH2I    = ThisGrid->GridData->AccessKDissH2I();
    int i, size=1;
    for (int dim=0; dim<rank; dim++)  size *= ArrDims[dim];
    for (i=0; i<size; i++)  phHI[i] = 0.0;
    if (RadiativeTransferHydrogenOnly == FALSE) {
      for (i=0; i<size; i++)  phHeI[i]  = 0.0;
      for (i=0; i<size; i++)  phHeII[i] = 0.0;
    }
    for (i=0; i<size; i++)  photogamma[i] = 0.0;
    if (MultiSpecies > 1) 
      for (i=0; i<size; i++)  dissH2I[i] = 0.0;
  }

  return SUCCESS;
}

#endif
