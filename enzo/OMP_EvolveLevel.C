/***********************************************************************
/
/  EVOLVE LEVEL FUNCTION
/
/  written by: Greg Bryan
/  date:       November, 1994
/  modified1:  February, 1995 by GB
/              Overhauled to make sure that all the subgrid's of a grid
/              advance with in lock step (i.e. with the same timestep and
/              in order).  This was done to allow a subgrid to get it's
/              boundary values from another subgrid (with the same parent).
/              Previously, a subgrid' BVs were always interpolated from its
/              parent.
/  modified2:  August, 1995 by GB
/                1) All grids on a level are processed at the same time
/                 (rather than all the subgrids of one parent).
/                2) C routines are called to loop over subgrids
/                 (so parallelizing C compilers can be used).
/                3) Subgrid timesteps are not constant over top grid step.
/              June, 1999 by GB -- Clean up somewhat
/
/  modified3:  August, 2001 by Alexei Kritsuk
/                Added 2nd call of PrepareDensityField() to compute
/                grav. potential (to be written with other baryon fields).
/  modified4:  January, 2004 by Alexei Kritsuk
/                Added support for RandomForcing
/  modified5:  February, 2006 by Daniel Reynolds
/                Added PotentialBdry to EvolveLevel and 
/                PrepareDensityField calls, so that it can be used
/                within computing isolating BCs for self-gravity.
/  modified6:  January, 2007 by Robert Harkness
/                Group and in-core i/o
/  modified7:  December, 2007 by Robert Harkness
/                Optional StaticSiblingList for root grid
/  modified8:  August, 2008 by Robert Harkness
/                Hybrid parallelism with OpenMP
/  modified9:  March, 2010 by Geoffrey So & Robert Harkness
/                Add emissivity field & rad hydro
/
/  PURPOSE:
/    This routine is the main grid evolution function.  It assumes that the
/    grids of level-1 have already been advanced by dt (passed
/    in the argument) and that their boundary values are properly set.
/    We then perform a complete update on all grids on level, including:
/       - for each grid: set the boundary values from parent/subgrids
/       - for each grid: get a list of its subgrids
/       - determine the timestep based on the minimum timestep for all grids
/       - subcycle over the grid timestep and for each grid:
/           - copy the fields to the old fields
/           - solve the hydro equations (and save fluxes around subgrid)
/           - set the boundary values from parent and/or other grids
/           - update time and check dt(min) for that grid against dt(cycle)
/           - call EvolveLevel(level+1)
/           - accumulate flux around this grid
/       - correct the solution on this grid based on subgrid solutions
/       - correct the solution on this grid based on improved subgrid fluxes
/
/    This routine essentially solves (completely) the grids of this level
/       and all finer levels and then corrects the solution of
/       grids on this level based on the improved subgrid solutions/fluxes.
/
/    Note: as a convenience, we store the current grid's fluxes (i.e. the
/          fluxes around the exterior of this grid) as the last entry in
/          the list of subgrids.
/
************************************************************************/

#ifdef RAD_HYDRO
#include "ImplicitProblemABC_preincludes.h"
#endif
 
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#ifdef USE_MPI
#include <mpi.h>
#endif
 
#include "performance.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "Hierarchy.h"
#include "TopGridData.h"
#include "LevelHierarchy.h"
#ifdef RAD_HYDRO
#include "ImplicitProblemABC.h"
#endif
 
/* function prototypes */
 
void DeleteFluxes(fluxes *Fluxes);
int  RebuildHierarchy(TopGridData *MetaData,
		      LevelHierarchyEntry *LevelArray[], int level);
int  ReportMemoryUsage(char *header = NULL);
int  UpdateParticlePositions(grid *Grid);
int  CheckEnergyConservation(HierarchyEntry *Grids[], int grid,
			     int NumberOfGrids, int level, float dt);
float CommunicationMinValue(float Value);
int GenerateGridArray(LevelHierarchyEntry *LevelArray[], int level,
		      HierarchyEntry **Grids[]);
 
#ifdef SIB3
int PrepareDensityField(LevelHierarchyEntry *LevelArray[],
			SiblingGridList SiblingList[],
			int level, TopGridData *MetaData, FLOAT When);
#else  // !SIB3
int PrepareDensityField(LevelHierarchyEntry *LevelArray[],
                        int level, TopGridData *MetaData, FLOAT When);
#endif  // end SIB3
 
#ifdef SIB2
int SetBoundaryConditions(HierarchyEntry *Grids[], int NumberOfGrids,
			  SiblingGridList SiblingList[],
			  int level, TopGridData *MetaData,
			  ExternalBoundary *Exterior, LevelHierarchyEntry * Level);
#else
int SetBoundaryConditions(HierarchyEntry *Grids[], int NumberOfGrids,
                          int level, TopGridData *MetaData,
                          ExternalBoundary *Exterior, LevelHierarchyEntry * Level);
#endif

#ifdef SAB
#ifdef SIB2
int SetAccelerationBoundary(HierarchyEntry *Grids[], int NumberOfGrids,
			    SiblingGridList SiblingList[],
			    int level, TopGridData *MetaData,
			    ExternalBoundary *Exterior,
			    LevelHierarchyEntry * Level,
			    int CycleNumber);
#else
int SetAccelerationBoundary(HierarchyEntry *Grids[], int NumberOfGrids,
			    int level, TopGridData *MetaData, 
			    ExternalBoundary *Exterior,
			    LevelHierarchyEntry * Level,
			    int CycleNumber);
#endif
#endif

#ifdef FLUX_FIX
int UpdateFromFinerGrids(int level, HierarchyEntry *Grids[], int NumberOfGrids,
			 int NumberOfSubgrids[],
			 fluxes **SubgridFluxesEstimate[],
			 LevelHierarchyEntry *SUBlingList[],
			 TopGridData *MetaData);
#else
int UpdateFromFinerGrids(int level, HierarchyEntry *Grids[], int NumberOfGrids,
			 int NumberOfSubgrids[],
			 fluxes **SubgridFluxesEstimate[]);
#endif
 
int CommunicationUpdateStarParticleCount(HierarchyEntry *Grids[],
					 TopGridData *MetaData,
					 int NumberOfGrids);
int RadiationFieldUpdate(LevelHierarchyEntry *LevelArray[], int level,
			 TopGridData *MetaData);
int WriteStreamData(HierarchyEntry *Grids[], int NumberOfGrids, 
		    TopGridData *MetaData, int CycleCount, int EndStep = FALSE);
int WriteMovieData(char *basename, int filenumber,
		   LevelHierarchyEntry *LevelArray[], TopGridData *MetaData,
		   FLOAT WriteTime);
int WriteTracerParticleData(char *basename, int filenumber,
		   LevelHierarchyEntry *LevelArray[], TopGridData *MetaData,
		   FLOAT WriteTime);

#ifdef USE_HDF5_GROUPS
int Group_WriteAllData(char *basename, int filenumber, HierarchyEntry *TopGrid,
		 TopGridData &MetaData, ExternalBoundary *Exterior,
#ifdef RAD_HYDRO
		 ImplicitProblemABC *ImplicitSolver,
#endif
		 FLOAT WriteTime = -1);
#else
int WriteAllData(char *basename, int filenumber, HierarchyEntry *TopGrid,
		 TopGridData &MetaData, ExternalBoundary *Exterior,
#ifdef RAD_HYDRO
		 ImplicitProblemABC *ImplicitSolver,
#endif
		 FLOAT WriteTime = -1);
#endif
 
int ComputeRandomForcingNormalization(LevelHierarchyEntry *LevelArray[],
                                      int level, TopGridData *MetaData,
                                      float * norm, float * pTopGridTimeStep);

int FastSiblingLocatorInitializeStaticChainingMesh(ChainingMeshStructure *Mesh, int Rank,
						   int TopGridDims[]); 
int FastSiblingLocatorInitialize(ChainingMeshStructure *Mesh, int Rank,
				 int TopGridDims[]);
int FastSiblingLocatorFinalize(ChainingMeshStructure *Mesh);
 
#ifdef FLUX_FIX
int CreateSUBlingList(TopGridData *MetaData,
		      HierarchyEntry *Grids[],
		      int NumberOfGrids,
		      LevelHierarchyEntry ***SUBlingList);
int DeleteSUBlingList(int NumberOfGrids,
		      LevelHierarchyEntry **SUBlingList);
#endif

void my_exit(int status);

#ifdef MEM_TRACE
Eint64 mused(void);
#endif



 
static int LevelCycleCount[MAX_DEPTH_OF_HIERARCHY];
double LevelWallTime[MAX_DEPTH_OF_HIERARCHY];
double LevelZoneCycleCount[MAX_DEPTH_OF_HIERARCHY];
double LevelZoneCycleCountPerProc[MAX_DEPTH_OF_HIERARCHY];
 
static float norm = 0.0;            //AK
static float TopGridTimeStep = 0.0; //AK

static int StaticSiblingListInitialized = 0;

#ifdef STATIC_SIBLING_LIST
static SiblingGridList StaticSiblingList[MAX_NUMBER_OF_SUBGRIDS];
static int StaticLevelZero = 1;
#else
static int StaticLevelZero = 0;
#endif

 


int EvolveLevel(TopGridData *MetaData, LevelHierarchyEntry *LevelArray[],
		int level, float dtLevelAbove, ExternalBoundary *Exterior
#ifdef RAD_HYDRO
		, ImplicitProblemABC *ImplicitSolver
#endif
		)
{

  int dbx = 0;
 
  FLOAT When;
  float dtThisLevelSoFar = 0.0, dtThisLevel, dtGrid, dtActual, dtLimit;
  int RefinementFactors[MAX_DIMENSION];
  int cycle = 0, counter = 0, grid1, subgrid;
  int Dummy;

  int max_thread_level = 30;
  int loop_threading;
  int min_grid_count = 1;

  HierarchyEntry *NextGrid;

  double t_entry, t_exit, t_call, t_return, t_rebuild1, t_rebuild2, t_acc;

#ifdef MEM_TRACE
    Eint64 MemInUse;
#endif

#ifdef USE_MPI
  t_entry = MPI_Wtime();
  t_acc = 0.0;
#endif




//  BEGIN

#ifdef MEM_TRACE
    MemInUse = mused();
    fprintf(memtracePtr, "Enter EL Level %8"ISYM"  %16"ISYM" \n", level, MemInUse);
    fflush(memtracePtr);
#endif

  /* Create an array (Grids) of all the grids. */

  typedef HierarchyEntry* HierarchyEntryPointer;
  HierarchyEntry **Grids;

  int NumberOfGrids = GenerateGridArray(LevelArray, level, &Grids);

  int *NumberOfSubgrids = new int[NumberOfGrids];
  fluxes ***SubgridFluxesEstimate = new fluxes **[NumberOfGrids];
  int *omp_error_flags = new int[NumberOfGrids];
  int *status = new int[NumberOfGrids];

#ifdef FLUX_FIX
  /* Create a SUBling list of the subgrids */
 
  LevelHierarchyEntry **SUBlingList;
#endif

  /* Initialize the chaining mesh used in the FastSiblingLocator. */

  if (dbx) fprintf(stderr, "EL: Initialize FSL \n"); 

  // If this is level 0 the SiblingList does not change and can be static

#ifdef STATIC_SIBLING_LIST
  if ( StaticLevelZero == 1 && level == 0 ) {

    if (!StaticSiblingListInitialized) {

      fprintf(stderr, "INITIALIZE Level 0 StaticSiblingList\n");

      ChainingMeshStructure StaticChainingMesh;

      FastSiblingLocatorInitializeStaticChainingMesh(&StaticChainingMesh, MetaData->TopGridRank,
						     MetaData->TopGridDims);

      for (grid1 = 0; grid1 < NumberOfGrids; grid1++)
        Grids[grid1]->GridData->FastSiblingLocatorAddGrid(&StaticChainingMesh);

      for (grid1 = 0; grid1 < NumberOfGrids; grid1++) {
        status[grid1] = Grids[grid1]->GridData->FastSiblingLocatorFindSiblings(
                              &StaticChainingMesh, &StaticSiblingList[grid1],
                              MetaData->LeftFaceBoundaryCondition,
                              MetaData->RightFaceBoundaryCondition);
      }

      /* Clean up the chaining mesh. */

      FastSiblingLocatorFinalize(&StaticChainingMesh);

      StaticSiblingListInitialized = 1;

    }

  } // if StaticLevelZero && level == 0
#endif

  SiblingGridList *SiblingList = new SiblingGridList[NumberOfGrids];

  ChainingMeshStructure ChainingMesh;

#ifdef STATIC_SIBLING_LIST
  if (StaticLevelZero == 1 && level == 0 ) {

    for (grid1 = 0; grid1 < NumberOfGrids; grid1++) {
      SiblingList[grid1].NumberOfSiblings = StaticSiblingList[grid1].NumberOfSiblings;
      SiblingList[grid1].GridList = StaticSiblingList[grid1].GridList;
    }

  }
#endif

  if (( StaticLevelZero == 1 && level != 0 ) || StaticLevelZero == 0 ) {

  FastSiblingLocatorInitialize(&ChainingMesh, MetaData->TopGridRank,
			       MetaData->TopGridDims);
 
  /* Add all the grids to the chaining mesh. */

  if (dbx) fprintf(stderr, "EL: FSL AddGrid entry \n");

  for (grid1 = 0; grid1 < NumberOfGrids; grid1++)
    Grids[grid1]->GridData->FastSiblingLocatorAddGrid(&ChainingMesh);

  if (dbx) fprintf(stderr, "EL: FSL AddGrid exit \n");
 
  /* For each grid, get a list of possible siblings from the chaining mesh. */

  for (grid1 = 0; grid1 < NumberOfGrids; grid1++) {
    status[grid1] = Grids[grid1]->GridData->FastSiblingLocatorFindSiblings(
                              &ChainingMesh, &SiblingList[grid1],
			      MetaData->LeftFaceBoundaryCondition,
			      MetaData->RightFaceBoundaryCondition);
  }

  /* Clean up the chaining mesh. */
 
  FastSiblingLocatorFinalize(&ChainingMesh);

  }


#ifdef EMISSIVITY
/* reset Emissivity array here before next step calculate the new values */
  if (StarMakerEmissivityField > 0) {

  /*
     clear the Emissivity of the level below, after the level below
     updated the current level (it's parent) and before the next
     timestep at the current level.
  */
    LevelHierarchyEntry *Temp;
    Temp = LevelArray[level+1];

    while (Temp != NULL) {
      /*
      printf("=======at level %"ISYM" clearing level %"ISYM"=======\n",
	     level, level+1);
      */
      Temp->GridData->ClearEmissivity();
      Temp = Temp->NextGridThisLevel;
    }
  }
#endif




  /* ================================================================== */
  /* For each grid: a) interpolate boundaries from its parent.
                    b) copy any overlapping zones.  */

#ifdef SIB2
  if (SetBoundaryConditions(Grids, NumberOfGrids, SiblingList,
			    level, MetaData, Exterior, LevelArray[level]) == FAIL)
    return FAIL;
#else
  if (SetBoundaryConditions(Grids, NumberOfGrids, level, MetaData,
                            Exterior, LevelArray[level]) == FAIL)
    return FAIL;
#endif




  /* Clear the boundary fluxes for all Grids (this will be accumulated over
     the subcycles below (i.e. during one current grid step) and used to by the
     current grid to correct the zones surrounding this subgrid (step #18). */
 
#pragma omp parallel private(grid1) shared(Grids, NumberOfGrids) default(none)
  {

#pragma omp for schedule(dynamic)
    for (grid1 = 0; grid1 < NumberOfGrids; grid1++) {
      Grids[grid1]->GridData->ClearBoundaryFluxes();
    }

  } // End omp parallel

  for (grid1 = 0; grid1 < NumberOfGrids; grid1++)
    Grids[grid1]->GridData->CountGridsOnLevel(level);


//loop_threading = level > 0 && level < max_thread_level+1;

  loop_threading = (level > 0) && (level_count[level] > min_grid_count);


  /* ================================================================== */
  /* Loop over grid timesteps until the elapsed time equals the timestep
     from the level above (or loop once for the top level). */

  dtActual = 0.0;

  while (dtThisLevelSoFar < dtLevelAbove) {
 
    /* Determine the timestep for this iteration of the loop. */
 
    if (level == 0) {
 
      /* For root level, use dtLevelAbove. */
 
      dtThisLevel      = dtLevelAbove;
      dtThisLevelSoFar = dtLevelAbove;
 
    } else {
 
      /* Compute the mininum timestep for all grids. */
 
      dtThisLevel = huge_number;
      for (grid1 = 0; grid1 < NumberOfGrids; grid1++) {
	dtGrid      = Grids[grid1]->GridData->ComputeTimeStep();
	dtThisLevel = min(dtThisLevel, dtGrid);
      }
      dtThisLevel = CommunicationMinValue(dtThisLevel);

      dtActual = dtThisLevel;

#ifdef USE_DT_LIMIT

//    dtLimit = LevelZeroDeltaT/(4.0)/POW(RefineBy,level);

      dtLimit = 0.5/(4.0)/POW(2.0,level);

      if ( dtActual < dtLimit ) {
        dtThisLevel = dtLimit;
      }

#endif
 
      /* Advance dtThisLevelSoFar (don't go over dtLevelAbove). */
 
      if (dtThisLevelSoFar+dtThisLevel*1.05 >= dtLevelAbove) {
	dtThisLevel      = dtLevelAbove - dtThisLevelSoFar;
	dtThisLevelSoFar = dtLevelAbove;
      }
      else
	dtThisLevelSoFar += dtThisLevel;
 
    }

    if (debug) printf("Level[%"ISYM"]: dt = %"GSYM"  %"GSYM"  (%"GSYM"/%"GSYM")\n", level, dtThisLevel, dtActual,
		      dtThisLevelSoFar, dtLevelAbove);
 
    /* Set all grid's timestep to this minimum dt. */

#pragma omp parallel private(grid1) shared(Grids, NumberOfGrids, dtThisLevel) default(none)
  {

#pragma omp for schedule(dynamic)
    for (grid1 = 0; grid1 < NumberOfGrids; grid1++) {
      Grids[grid1]->GridData->SetTimeStep(dtThisLevel);
    }

  } // End omp parallel


    /* For each grid, compute the number of it's subgrids. */
 
//PRAG omp parallel private(grid1, counter, NextGrid) shared(Grids, NumberOfGrids, NumberOfSubgrids) default(none)
//   {
//PRAG omp for schedule(dynamic)
    for (grid1 = 0; grid1 < NumberOfGrids; grid1++) {

      NextGrid = Grids[grid1]->NextGridNextLevel;
      counter = 0;

      while (NextGrid != NULL) {
	NextGrid = NextGrid->NextGridThisLevel;
	if (++counter > MAX_NUMBER_OF_SUBGRIDS) {
	  fprintf(stderr, "More subgrids than MAX_NUMBER_OF_SUBGRIDS.\n");
	  return FAIL;
	}
      }

      NumberOfSubgrids[grid1] = counter + 1;

    }

//PRAG  }  // End omp parallel
 

    /* For each grid, create the subgrid list. */
 
//PRAG? omp parallel private(grid1) shared(SubgridFluxesEstimate) default(none)
//  {
//PRAG? omp for schedule(dynamic)
    for (grid1 = 0; grid1 < NumberOfGrids; grid1++) {

      /* Allocate the subgrid fluxes for this grid. */

      SubgridFluxesEstimate[grid1] = new fluxes *[NumberOfSubgrids[grid1]];

    }
//  }

#pragma omp parallel private(grid1, subgrid) shared(NumberOfGrids, NumberOfSubgrids, SubgridFluxesEstimate) default(none)
  {

#pragma omp for schedule(dynamic)
    for (grid1 = 0; grid1 < NumberOfGrids; grid1++) {

      for (subgrid = 0; subgrid < NumberOfSubgrids[grid1]; subgrid++)
        SubgridFluxesEstimate[grid1][subgrid] = NULL;

    }

  } // End omp parallel

//PRAG? omp parallel private(grid1, counter, NextGrid, RefinementFactors) shared(MyProcessorNumber, Grids, NumberOfGrids, SubgridFluxesEstimate) default(none)
//  {
//PRAG? omp for schedule(dynamic)
    for (grid1 = 0; grid1 < NumberOfGrids; grid1++) {
 
      /* Collect the flux data and store it in the newly minted fluxes.
	 Or rather that's what we should do.  Instead, we create fluxes one
	 by one in this awkward array of pointers to pointers.  This should be
	 changed so that all the routines take arrays of flux rather than
	 arrays of pointers to flux.  Dumb. */
 
      counter = 0;

      if (MyProcessorNumber == Grids[grid1]->GridData->ReturnProcessorNumber()) {
 
	NextGrid = Grids[grid1]->NextGridNextLevel;

	while (NextGrid != NULL) {
	  SubgridFluxesEstimate[grid1][counter] = new fluxes;
	  Grids[grid1]->GridData->ComputeRefinementFactors(NextGrid->GridData, RefinementFactors);
	  NextGrid->GridData->ReturnFluxDims(*(SubgridFluxesEstimate[grid1][counter++]), RefinementFactors);
	  NextGrid = NextGrid->NextGridThisLevel;
	}
 
	/* Add the external boundary of this subgrid to the subgrid list. This
	   makes it easy to keep adding up the fluxes of this grid, but we must
	   keep in mind that the last subgrid should be ignored elsewhere. */
 
	SubgridFluxesEstimate[grid1][counter] = new fluxes;
	Grids[grid1]->GridData->ComputeRefinementFactors(Grids[grid1]->GridData, RefinementFactors);
	Grids[grid1]->GridData->ReturnFluxDims(*(SubgridFluxesEstimate[grid1][counter]), RefinementFactors);

      }
 
    } // end loop over grids (create Subgrid list)
//  }
 

    /* ------------------------------------------------------- */
    /* Prepare the density field (including particle density). */
 
//  fprintf(stderr, "%"ISYM": EvolveLevel: Enter PrepareDensityField\n", MyProcessorNumber);
 
    When = 0.5;

#ifdef RAD_HYDRO
    if (RadiationHydrodynamics < 2) {
#endif

#ifdef SIB3
    if (SelfGravity)
      if (PrepareDensityField(LevelArray, SiblingList, level, MetaData, When) == FAIL) {
	fprintf(stderr, "Error in PrepareDensityField.\n");
	return FAIL;
      }
#else   // !SIB3
    if (SelfGravity)
      if (PrepareDensityField(LevelArray, level, MetaData, When) == FAIL) {
        fprintf(stderr, "Error in PrepareDensityField.\n");
        return FAIL;
      }
#endif  // end SIB3
 
 
//  fprintf(stderr, "%"ISYM": EvolveLevel: Exit PrepareDensityField\n", MyProcessorNumber);
 
    /* Prepare normalization for random forcing. Involves top grid only. */
 
    if (RandomForcing && MetaData->CycleNumber > 0 && level == 0) {
      if ( ComputeRandomForcingNormalization(LevelArray, 0, MetaData, &norm, &TopGridTimeStep) == FAIL ) {
        fprintf(stderr, "Error in ComputeRandomForcingNormalization.\n");
        return FAIL;
      }
    }

#ifdef RAD_HYDRO
    }
#endif
 

    /* ------------------------------------------------------- */
    /* Evolve all grids by timestep dtThisLevel. */
 
    for (grid1 = 0; grid1 < NumberOfGrids; grid1++) {
 
      /* Call analysis routines. */
 
      if (ProblemType == 24)
	Grids[grid1]->GridData->SphericalInfallGetProfile(level, 1);

      if (ProblemType == 30)
	Grids[grid1]->GridData->AnalyzeTrackPeaks(level, 0);

      if (ProblemType == 27)
	if (Grids[grid1]->GridData->ReturnProcessorNumber() == MyProcessorNumber) {

	  float AM[3], MeanVelocity[3], DMVelocity[3];
	  FLOAT Center[] = {0,0,0}, CenterOfMass[3], DMCofM[3];

	  Grids[grid1]->GridData->CalculateAngularMomentum(Center, AM,
			   MeanVelocity, DMVelocity, CenterOfMass, DMCofM);

	  fprintf(stdout, "level = %"ISYM" %"ISYM" %"ISYM"  Vel %"FSYM" %"FSYM" %"FSYM"  DMVel %"FSYM" %"FSYM" %"FSYM"  CofM %"PSYM" %"PSYM" %"PSYM"  DMCofM %"FSYM" %"FSYM" %"FSYM"\n",
		level, LevelCycleCount[level], grid1, MeanVelocity[0],
		MeanVelocity[1], MeanVelocity[2],
		DMVelocity[0], DMVelocity[1], DMVelocity[2],
		-CenterOfMass[0], -CenterOfMass[1], -CenterOfMass[2],
		DMCofM[0], DMCofM[1], DMCofM[2]);
	}
 
    }


#ifdef RAD_HYDRO
    if (RadiationHydrodynamics < 2) {
#endif

#pragma omp parallel private(grid1, Dummy) shared(NumberOfGrids, Grids, level, MaximumGravityRefinementLevel, SelfGravity, status) default(none)
  {
#pragma omp for schedule(dynamic)
    for (grid1 = 0; grid1 < NumberOfGrids; grid1++) {

      /* Gravity: compute acceleration field for grid and particles. */
 
      if (SelfGravity) {

	if (level <= MaximumGravityRefinementLevel) {
 
	  /* Compute the potential. */
 
	  if (level > 0) {
	    status[grid1] = Grids[grid1]->GridData->SolveForPotential(Dummy, level);
	  }

	}

	  /* otherwise, interpolate potential from coarser grid, which is now done in PrepareDensity. */

      } // end: if (SelfGravity)

    } // End of loop over grids

  } // End parallel region

#ifdef RAD_HYDRO
    }
#endif


#ifdef RAD_HYDRO
    if (RadiationHydrodynamics < 2) {
#endif

    for (grid1 = 0; grid1 < NumberOfGrids; grid1++) {

      if (SelfGravity) {

	if (level <= MaximumGravityRefinementLevel) {

	  if (Grids[grid1]->GridData->ComputeAccelerations(level) == FAIL) {
	    fprintf(stderr, "Error in grid->ComputeAccelerations.\n");
	    return FAIL;
	  }

	}

      } // end: if (SelfGravity)
 
      /* Gravity: compute field due to preset sources. */
 
      if (UniformGravity || PointSourceGravity) {
	if (Grids[grid1]->GridData->ComputeAccelerationFieldExternal() == FAIL) {
	  fprintf(stderr,"Error in grid->ComputeAccelerationFieldExternal.\n");
	  return FAIL;
	}
      }
 
      /* Check for energy conservation. */

/*
      if (ComputePotential)
	if (CheckEnergyConservation(Grids, grid, NumberOfGrids, level, dtThisLevel) == FAIL) {
	  fprintf(stderr, "Error in CheckEnergyConservation.\n");
	  return FAIL;
	}
*/

    } // End of loop over grids

#ifdef RAD_HYDRO
    }
#endif



    //This ensures that all subgrids agree in the boundary.
    //Not a big deal for hydro, but essential for DivB = 0 in MHD runs.
    //Only called on level > 0 because the root grid is dealt with differently than SG's.

#ifdef RAD_HYDRO
    if (RadiationHydrodynamics < 2) {
#endif

#ifdef SAB
    if ( (SelfGravity || UniformGravity || PointSourceGravity) && level > 0) {

#ifdef SIB2
      if( SetAccelerationBoundary(Grids, NumberOfGrids,
				  SiblingList,
				  level, MetaData,
				  Exterior, LevelArray[level], LevelCycleCount[level]) == FAIL ) {
	fprintf(stderr,"Error with AccelerationBoundary.\n");
	return FAIL;
      }
#else
      if( SetAccelerationBoundary(Grids, NumberOfGrids,
				  level, MetaData,
				  Exterior, LevelArray[level], LevelCycleCount[level]) == FAIL ) {
	fprintf(stderr,"Error with AccelerationBoundary.\n");
	return FAIL;
      }
#endif

    }
#endif

#ifdef RAD_HYDRO
    }
#endif



//  OpenMP Loop over grids

#ifdef RAD_HYDRO
#pragma omp parallel private(grid1) shared(level, LevelCycleCount, NumberOfGrids, NumberOfSubgrids, Grids, MetaData, RandomForcing, TopGridTimeStep, norm, SubgridFluxesEstimate, MultiSpecies, RadiativeCooling, StarParticleCreation, StarParticleFeedback, SelfGravity, UniformGravity, PointSourceGravity, MaximumGravityRefinementLevel, MaximumRefinementLevel, LevelArray, ComovingCoordinates, RadiationHydrodynamics, RadiativeTransfer, ImplicitSolver, dtLevelAbove, status) default(none) if(loop_threading)
  {
#else
#pragma omp parallel private(grid1) shared(level, LevelCycleCount, NumberOfGrids, NumberOfSubgrids, Grids, MetaData, RandomForcing, TopGridTimeStep, norm, SubgridFluxesEstimate, MultiSpecies, RadiativeCooling, StarParticleCreation, StarParticleFeedback, SelfGravity, UniformGravity, PointSourceGravity, MaximumGravityRefinementLevel, MaximumRefinementLevel, LevelArray, ComovingCoordinates, dtLevelAbove, status) default(none) if(loop_threading)
  {
#endif

#pragma omp for schedule(dynamic)
    for (grid1 = 0; grid1 < NumberOfGrids; grid1++) {

      /* Copy current fields (with their boundaries) to the old fields
	  in preparation for the new step. */
 
      status[grid1] = Grids[grid1]->GridData->CopyBaryonFieldToOldBaryonField();

    }


#ifdef RAD_HYDRO
    if (RadiationHydrodynamics < 2) {
#endif

#pragma omp for schedule(dynamic)
    for (grid1 = 0; grid1 < NumberOfGrids; grid1++) {

      /* Add RandomForcing fields to velocities after the copying of current
         fields to old. I also update the total energy accordingly here.
         It makes no sense to force on the very first time step. */
 
      if (RandomForcing && MetaData->CycleNumber > 0) //AK
        status[grid1] = Grids[grid1]->GridData->AddRandomForcing(&norm, TopGridTimeStep);
 
    }

#ifdef RAD_HYDRO
    }
#endif




#ifdef RAD_HYDRO
    if (RadiationHydrodynamics < 2) {
#endif

#pragma omp for schedule(dynamic)
    for (grid1 = 0; grid1 < NumberOfGrids; grid1++) {

      /* Call hydro solver and save fluxes around subgrids. */
 
//      fprintf(stderr, "%"ISYM": Calling Hydro\n", MyProcessorNumber);
 
      status[grid1] =Grids[grid1]->GridData->SolveHydroEquations(LevelCycleCount[level],
       	 NumberOfSubgrids[grid1], SubgridFluxesEstimate[grid1], level);
 
//      fprintf(stderr, "%"ISYM": Called Hydro\n", MyProcessorNumber);

    }

#ifdef RAD_HYDRO
    }
#endif




//-----------------------------------------------------------------------------------------


#ifdef RAD_HYDRO
    if ((RadiationHydrodynamics == 0) || (level > 0) || (RadiativeTransfer > 0)) {
#endif

#ifdef RATE_AND_COOL

    if (MultiSpecies && RadiativeCooling) {

#pragma omp for schedule(dynamic)
    for (grid1 = 0; grid1 < NumberOfGrids; grid1++) {

      /* Solve the cooling and species rate equations. */
 
//      fprintf(stderr, "%"ISYM": Calling SolveCoolAndRateEquations\n", MyProcessorNumber);
 
	status[grid1] = Grids[grid1]->GridData->SolveRateAndCoolEquations();
 
//      fprintf(stderr, "%"ISYM": Called SolveCoolAndRateEquations\n", MyProcessorNumber);

    }

    } else {

#endif

//-----------------------------------------------------------------------------------------

#pragma omp for schedule(dynamic)
    for (grid1 = 0; grid1 < NumberOfGrids; grid1++) {

//      fprintf(stderr, "%"ISYM": Calling MultiSpecies\n", MyProcessorNumber);
 
	if (MultiSpecies)
          status[grid1] = Grids[grid1]->GridData->SolveRateEquations();
 
//      fprintf(stderr, "%"ISYM": Called MultiSpecies\n", MyProcessorNumber);

    }


#pragma omp for schedule(dynamic)
    for (grid1 = 0; grid1 < NumberOfGrids; grid1++) {
 
	/* Include radiative cooling/heating. */
 
//      fprintf(stderr, "%"ISYM": Calling RadiativeCooling\n", MyProcessorNumber);
 
	if (RadiativeCooling)
          status[grid1] = Grids[grid1]->GridData->SolveRadiativeCooling();
 
//      fprintf(stderr, "%"ISYM": Called RadiativeCooling\n", MyProcessorNumber);

    }

#ifdef RATE_AND_COOL
    }
#endif

//-----------------------------------------------------------------------------------------

#ifdef RAD_HYDRO
      }  /* end if (RadiationHydrodynamics == 0 || level > 0 || RadiativeTransfer > 0) */
#endif



#ifdef RAD_HYDRO
      if (RadiationHydrodynamics < 2) {
#endif

#pragma omp for schedule(dynamic)
    for (grid1 = 0; grid1 < NumberOfGrids; grid1++) {
 
      /* Update particle positions (if present). */
 
//      fprintf(stderr, "%"ISYM": Calling UpdatePP\n", MyProcessorNumber);
 
      status[grid1] = UpdateParticlePositions(Grids[grid1]->GridData);
 
//      fprintf(stderr, "%"ISYM": Called UpdatePP\n", MyProcessorNumber);

    }

#ifdef RAD_HYDRO
      }
#endif


#ifdef RAD_HYDRO
      if (RadiationHydrodynamics < 2) {
#endif

#pragma omp for schedule(dynamic)
    for (grid1 = 0; grid1 < NumberOfGrids; grid1++) {
 
      /* Include 'star' particle creation and feedback.
         (first, set the under_subgrid field). */
 
      if (StarParticleCreation || StarParticleFeedback) {

	status[grid1] = Grids[grid1]->GridData->ZeroSolutionUnderSubgrid(NULL,
						 ZERO_UNDER_SUBGRID_FIELD);

	LevelHierarchyEntry *Temp2 = LevelArray[level+1];

	while (Temp2 != NULL) {
	  Grids[grid1]->GridData->ZeroSolutionUnderSubgrid(Temp2->GridData,
					 ZERO_UNDER_SUBGRID_FIELD);
	  Temp2 = Temp2->NextGridThisLevel;
	}

      }

      if (StarParticleCreation || StarParticleFeedback) {

#ifdef EMISSIVITY
	status[grid1] = Grids[grid1]->GridData->StarParticleHandler(level, dtLevelAbove);
#else
	status[grid1] = Grids[grid1]->GridData->StarParticleHandler(level);
#endif

      }
 
    }

#ifdef RAD_HYDRO
      }
#endif

      /* Coupled Rad-Hydro solver. */

#ifdef RAD_HYDRO
      if ((RadiationHydrodynamics > 0) && (level == 0))  {

        for (grid1 = 0; grid1 < NumberOfGrids; grid1++) {
          status[grid1] = ImplicitSolver->Evolve(Grids[grid1], Grids[grid1]->GridData->ReturnTimeStep());
	}

      }
#endif



#ifdef RAD_HYDRO
      if (RadiationHydrodynamics < 2) {
#endif

#pragma omp for schedule(dynamic)
    for (grid1 = 0; grid1 < NumberOfGrids; grid1++) {

      /* Gravity: clean up AccelerationField. */

      // David Collins removes this for MHD Amr
      if (SelfGravity || UniformGravity || PointSourceGravity) {
	if (level != MaximumGravityRefinementLevel ||
	    MaximumGravityRefinementLevel == MaximumRefinementLevel)
	  Grids[grid1]->GridData->DeleteAccelerationField();
	Grids[grid1]->GridData->DeleteParticleAcceleration();
      }
 
    }

#ifdef RAD_HYDRO
      }
#endif




#pragma omp for schedule(dynamic)
    for (grid1 = 0; grid1 < NumberOfGrids; grid1++) {

      /* Update current problem time of this subgrid. */
 
      Grids[grid1]->GridData->SetTimeNextTimestep();
 
    }

#pragma omp for schedule(dynamic)
    for (grid1 = 0; grid1 < NumberOfGrids; grid1++) {

      /* If using comoving co-ordinates, do the expansion terms now. */
 
      if (ComovingCoordinates)
	status[grid1] = Grids[grid1]->GridData->ComovingExpansionTerms();
 
    }  // end loop over grids

  } // End parallel region

// End OpenMP Parallel Region



 
    /* For each grid: a) interpolate boundaries from the parent grid.
                      b) copy any overlapping zones from siblings. */
 
#ifdef SIB2
    if (SetBoundaryConditions(Grids, NumberOfGrids, SiblingList,
			      level, MetaData, Exterior, LevelArray[level]) == FAIL)
      return FAIL;
#else
    if (SetBoundaryConditions(Grids, NumberOfGrids, level, MetaData,
                              Exterior, LevelArray[level]) == FAIL)
      return FAIL;
#endif




    /* Update the star particle counters. */
 
    if (StarParticleCreation) {
      if (CommunicationUpdateStarParticleCount(Grids, MetaData, NumberOfGrids) == FAIL)
	return FAIL;
    }
 
    /* Check for movie output (only check if this is bottom of hierarchy). */
 
    if (LevelArray[level+1] == NULL) {
      if (LevelArray[level]->GridData->ReturnTime() >=
	  MetaData->TimeLastMovieDump + MetaData->dtMovieDump &&
	  MetaData->dtMovieDump > 0.0) {
	MetaData->TimeLastMovieDump += MetaData->dtMovieDump;
	if (WriteMovieData(MetaData->MovieDumpName,
			  MetaData->MovieDumpNumber++, LevelArray, MetaData,
			  LevelArray[level]->GridData->ReturnTime()) == FAIL) {
	  fprintf(stderr, "Error in WriteMovieData.\n");
	  return FAIL;
	}
      }
    }
 
    /* Check for tracer particle output (only if this bottom of hierarchy). */
 
    if (LevelArray[level+1] == NULL) {
      if (LevelArray[level]->GridData->ReturnTime() >=
	  MetaData->TimeLastTracerParticleDump +
	  MetaData->dtTracerParticleDump &&
	  MetaData->dtTracerParticleDump > 0.0) {
	MetaData->TimeLastTracerParticleDump += MetaData->dtTracerParticleDump;
	if (WriteTracerParticleData(MetaData->TracerParticleDumpName,
				    MetaData->TracerParticleDumpNumber++,
				    LevelArray, MetaData,
			  LevelArray[level]->GridData->ReturnTime()) == FAIL) {
	  fprintf(stderr, "Error in WriteTracerParticleData.\n");
	  return FAIL;
	}
      }
    }
 



    /* If cosmology, then compute grav. potential for output if needed. */

#ifdef RAD_HYDRO
    if (RadiationHydrodynamics < 2) {
#endif
 
    if (ComovingCoordinates && SelfGravity && WritePotential) {

      CopyGravPotential = TRUE;
      When = 0.0;
 
#ifdef SIB3
      if (PrepareDensityField(LevelArray, SiblingList, level, MetaData, When) == FAIL) {
        fprintf(stderr, "Error in PrepareDensityField.\n");
        return FAIL;
      }
#else   // !SIB3
      if (PrepareDensityField(LevelArray, level, MetaData, When) == FAIL) {
        fprintf(stderr, "Error in PrepareDensityField.\n");
        return FAIL;
      }
#endif  // end SIB3
 
      CopyGravPotential = FALSE;

#pragma omp parallel private(grid1, Dummy) shared(NumberOfGrids, Grids, level, MaximumGravityRefinementLevel, status) default(none)
  {
#pragma omp for schedule(dynamic)
      for (grid1 = 0; grid1 < NumberOfGrids; grid1++) {

        if (level <= MaximumGravityRefinementLevel) {
 
          /* Compute the potential. */
 
          if (level > 0) {
            status[grid1] = Grids[grid1]->GridData->SolveForPotential(Dummy, level);
          }

        }

      } //  end loop over grids

  } // end parallel region

//PRAG? omp parallel private(grid1) shared(NumberOfGrids, Grids, level, MaximumGravityRefinementLevel) default(none)
//  {
//PRAG? omp for schedule(dynamic)
      for (grid1 = 0; grid1 < NumberOfGrids; grid1++) {
        if (level <= MaximumGravityRefinementLevel) {
          Grids[grid1]->GridData->CopyPotentialToBaryonField();
        }
      }
//  }

    } // if WritePotential
 
#ifdef RAD_HYDRO
    }
#endif



    /* Check for new level output (only if this is bottom of hierarchy). */
 
    if (MetaData->OutputFirstTimeAtLevel > 0 &&
	level >= MetaData->OutputFirstTimeAtLevel &&
	LevelArray[level+1] == NULL) {

      MetaData->OutputFirstTimeAtLevel = level+1;

      LevelHierarchyEntry *Temp2 = LevelArray[0];

      while (Temp2->NextGridThisLevel != NULL)
	Temp2 = Temp2->NextGridThisLevel; /* ugh: find last in linked list */

#ifdef USE_HDF5_GROUPS
      if (Group_WriteAllData(MetaData->DataDumpName, MetaData->DataDumpNumber++,
		       Temp2->GridHierarchyEntry, *MetaData, Exterior,
#ifdef RAD_HYDRO
		       ImplicitSolver,
#endif
		       LevelArray[level]->GridData->ReturnTime()) == FAIL) {
	fprintf(stderr, "Error in Group_WriteAllData.\n");
	return FAIL;
      }
#else
      if (WriteAllData(MetaData->DataDumpName, MetaData->DataDumpNumber++,
		       Temp2->GridHierarchyEntry, *MetaData, Exterior, 
#ifdef RAD_HYDRO
		       ImplicitSolver,
#endif
		       LevelArray[level]->GridData->ReturnTime()) == FAIL) {
	fprintf(stderr, "Error in WriteAllData.\n");
	return FAIL;
      }
#endif

    } // OutputFirstTimeAtLevel
 
    /* Check for stop (unpleasant to exit from here, but...). */
 
    if (MetaData->StopFirstTimeAtLevel > 0 &&
	level >= MetaData->StopFirstTimeAtLevel &&
	LevelArray[level+1] == NULL) {

      // Write movie data in all grids if necessary

      if (MovieSkipTimestep != INT_UNDEFINED)
	for (int mlevel = 0; mlevel < MAX_DEPTH_OF_HIERARCHY; mlevel++) {
	  if (LevelArray[mlevel] == NULL) break;
	  delete [] Grids;
	  NumberOfGrids = GenerateGridArray(LevelArray, mlevel, &Grids);
	  if (WriteStreamData(Grids, NumberOfGrids, MetaData,
			      LevelCycleCount[mlevel], TRUE) == FAIL) {
	    fprintf(stderr, "Error in WriteStreamData.\n");
	    return FAIL;
	  }
	}

      fprintf(stderr, "Stopping due to request on level %"ISYM"\n", level);
      my_exit(EXIT_SUCCESS);

    } // StopFirstTimeAtLevel



 
    /* For each grid, delete the GravitatingMassFieldParticles. */

#ifdef RAD_HYDRO
    if (RadiationHydrodynamics < 2) {
#endif
 
#pragma omp parallel private(grid1) shared(Grids, NumberOfGrids) default(none)
  {

#pragma omp for schedule(dynamic)
    for (grid1 = 0; grid1 < NumberOfGrids; grid1++) {
      Grids[grid1]->GridData->DeleteGravitatingMassFieldParticles();
    }

  } // End omp parallel
 
#ifdef RAD_HYDRO
    }
#endif



    /* ----------------------------------------- */
    /* Evolve the next level down (recursively). */
 
    if (dbx) fprintf(stderr, "EL Level %"ISYM" going to Level %"ISYM"\n", level, level+1);

#ifdef MEM_TRACE
    MemInUse = mused();
    fprintf(memtracePtr, "  EL Level %"ISYM" to Level %"ISYM" %16"ISYM" \n", level, level+1, MemInUse);
    fflush(memtracePtr);
#endif

#ifdef USE_MPI
  t_call = MPI_Wtime();
#endif

    if (LevelArray[level+1] != NULL) {
      if (EvolveLevel(MetaData, LevelArray, level+1, dtThisLevel, Exterior
#ifdef RAD_HYDRO
		      , ImplicitSolver
#endif
		      ) == FAIL) {
	fprintf(stderr, "Error in EvolveLevel (%"ISYM").\n", level);
	return FAIL;
      }
    }

#ifdef USE_MPI
  t_return = MPI_Wtime();
#endif

#ifdef MEM_TRACE
    MemInUse = mused();
    fprintf(memtracePtr, "  EL Level %"ISYM" from Level %"ISYM" %16"ISYM" \n", level, level+1, MemInUse);
    fflush(memtracePtr);
#endif

    if (dbx) fprintf(stderr, "EL Level %"ISYM" returns from Level %"ISYM"\n", level, level+1);




    // Streaming movie output (only run if everything is evolved)

    if (MovieSkipTimestep != INT_UNDEFINED) {
      if (WriteStreamData(Grids, NumberOfGrids, MetaData, 
			  LevelCycleCount[level]) == FAIL) {
	fprintf(stderr, "Error in WriteStreamData.\n");
	return FAIL;
      }
    }




    /* ------------------------------------------------------- */
    /* For each grid,
     * (a) project the subgrid's solution into this grid (step #18)
     * (b) correct for the difference between this grid's fluxes and the
     *     subgrid's fluxes. (step #19)
     */
 
#ifdef FLUX_FIX

    SUBlingList = new LevelHierarchyEntry*[NumberOfGrids];

    for(int list=0; list < NumberOfGrids; list++)
      SUBlingList[list] = NULL;

 
    if (FluxCorrection) {

      /* Fill in the SUBling list */

      if (dbx) fprintf(stderr, "EL: CSL entry \n");
      if (CreateSUBlingList(MetaData, Grids,
                              NumberOfGrids, &SUBlingList) == FAIL) {
        fprintf(stderr, "Error in CreateSUBlingList.\n");
        return FAIL;
      }
      if (dbx) fprintf(stderr, "EL: CSL exit \n");

    }

#endif
 
#ifdef FLUX_FIX
    if (UpdateFromFinerGrids(level, Grids, NumberOfGrids, NumberOfSubgrids,
			     SubgridFluxesEstimate,
			     SUBlingList,
			     MetaData) == FAIL)
      return FAIL;
#else
    if (UpdateFromFinerGrids(level, Grids, NumberOfGrids, NumberOfSubgrids,
			     SubgridFluxesEstimate) == FAIL)
      return FAIL;
#endif

    if (dbx) fprintf(stderr, "OK after UpdateFromFinerGrids \n");

#ifdef FLUX_FIX
    if ( FluxCorrection ) {
      /* Clean up SUBlings */
      if (DeleteSUBlingList( NumberOfGrids, SUBlingList ) == FAIL) {
        fprintf(stderr, "Error in DeleteSUBlingList.\n");
        return FAIL;
      }
    }
#endif

    if (dbx) fprintf(stderr, "OK after DeleteSUBlingList \n");



 
    /* ------------------------------------------------------- */
    /* Add the saved fluxes (in the last subsubgrid entry) to the exterior
       fluxes for this subgrid.
       (Note: this must be done after CorrectForRefinedFluxes). */
 
    // NOTE: This is a reduction operation!

#ifdef RAD_HYDRO
    if (RadiationHydrodynamics < 2) {
#endif

#pragma omp parallel private(grid1) shared(MyProcessorNumber, NumberOfGrids, NumberOfSubgrids, Grids, SubgridFluxesEstimate, FluxCorrection, status) default(none)
  {
#pragma omp for schedule(dynamic)
    for (grid1 = 0; grid1 < NumberOfGrids; grid1++) {

      if (MyProcessorNumber == Grids[grid1]->GridData->ReturnProcessorNumber()) {
 
      if (FluxCorrection) {
	status[grid1] = Grids[grid1]->GridData->AddToBoundaryFluxes
	    (SubgridFluxesEstimate[grid1][NumberOfSubgrids[grid1] - 1]);
      }

      }

    } // end of loop over grids

  } // End parallel region


//PRAG? omp parallel private(grid1, subgrid) shared(MyProcessorNumber, NumberOfGrids, NumberOfSubgrids, SubgridFluxesEstimate) default(none)
//  {
//PRAG? omp for schedule(dynamic)
    for (grid1 = 0; grid1 < NumberOfGrids; grid1++) {

      if (MyProcessorNumber == Grids[grid1]->GridData->ReturnProcessorNumber()) {
 
      /* Delete fluxes pointed to by SubgridFluxesEstimate[subgrid]. */
 
      for (subgrid = 0; subgrid < NumberOfSubgrids[grid1]; subgrid++) {
	DeleteFluxes(SubgridFluxesEstimate[grid1][subgrid]);
	delete       SubgridFluxesEstimate[grid1][subgrid];
      }

      delete [] SubgridFluxesEstimate[grid1];

      }
 
    } // end of loop over grids

//  } 

#ifdef RAD_HYDRO
    }
#endif



    /* Recompute radiation field, if requested. */

#ifdef RAD_HYDRO
    if (RadiationHydrodynamics == 0) {
#endif
 
    if (RadiationFieldType >= 10 && RadiationFieldType <= 11 &&
	level <= RadiationFieldLevelRecompute)
      if (RadiationFieldUpdate(LevelArray, level, MetaData) == FAIL) {
	fprintf(stderr, "Error in RecomputeRadiationField.\n");
	return FAIL;
      }
 
#ifdef RAD_HYDRO
    }
#endif



    /* Rebuild the Grids on the next level down.
       Don't bother on the last cycle, as we'll rebuild this grid soon. */
 
#ifdef MEM_TRACE
    MemInUse = mused();
    fprintf(memtracePtr, "EL rebuild on level %8"ISYM"  %16"ISYM" \n", level, MemInUse);
    fflush(memtracePtr);
#endif

#ifdef USE_MPI
    t_rebuild1 = MPI_Wtime();
#endif

    if (dtThisLevelSoFar < dtLevelAbove) {
      if (RebuildHierarchy(MetaData, LevelArray, level) == FAIL) {
	fprintf(stderr, "Error in RebuildHierarchy.\n");
	return FAIL;
      }
    }

#ifdef USE_MPI
    t_rebuild2 = MPI_Wtime();
#endif

    t_acc += max((t_rebuild2 - t_rebuild1), 0.0);

#ifdef MEM_TRACE
    MemInUse = mused();
    fprintf(memtracePtr, "EL rebuild on level %8"ISYM"  %16"ISYM" \n", level, MemInUse);
    fflush(memtracePtr);
#endif
 
    /* Count up number of grids on this level. */
 
    int GridMemory, NumberOfCells, CellsTotal, Particles;
    float AxialRatio, GridVolume;

    for (grid1 = 0; grid1 < NumberOfGrids; grid1++) {
      Grids[grid1]->GridData->CollectGridInformation
        (GridMemory, GridVolume, NumberOfCells, AxialRatio, CellsTotal, Particles);
      LevelZoneCycleCount[level] += NumberOfCells;
      if (MyProcessorNumber == Grids[grid1]->GridData->ReturnProcessorNumber())
	LevelZoneCycleCountPerProc[level] += NumberOfCells;
    }
 
    cycle++;
    LevelCycleCount[level]++;
 
  } // end of while loop over subcycles



 
  if (debug)
    fprintf(stderr, "EvolveLevel[%"ISYM"]: NumberOfSubCycles = %"ISYM" (%"ISYM" total)\n",
            level, cycle, LevelCycleCount[level]);
 
  /* Clean up. */
 
  delete [] NumberOfSubgrids;
  delete [] Grids;
  delete [] SubgridFluxesEstimate;
 
  /* Clean up the sibling list. */

  if (( StaticLevelZero == 1 && level != 0 ) || StaticLevelZero == 0 ) {
    for (grid1 = 0; grid1 < NumberOfGrids; grid1++)
      delete [] SiblingList[grid1].GridList;
    delete [] SiblingList;
  }

#ifdef MEM_TRACE
    MemInUse = mused();
    fprintf(memtracePtr, "Exit EL level %8"ISYM"  %16"ISYM" \n", level, MemInUse);
    fflush(memtracePtr);
#endif

#ifdef USE_MPI
  t_exit = MPI_Wtime();
#endif

  rebuild_timer[level] += t_acc;
  level_timer[level] += (max((t_exit-t_entry), 0.0) - max((t_return-t_call), 0.0));

  if (dbx) fprintf(stderr, "Return from EL Level %"ISYM"\n", level);
 
  return SUCCESS;
 
}
