/***********************************************************************
/
/  EVOLVE LEVEL ROUTINES (CALLED BY EVOLVE LEVEL)
/
/  written by: Greg Bryan
/  date:       June, 1999
/  modifiedN:  Robert Harkness
/  date:       29th October, 2008
/
/  PURPOSE:  This is a collection of routines called by EvolveLevel.
/            These have been optimized for enhanced message passing
/            performance by performing two passes -- one which generates
/            sends and the second which receives them.
/
/  modified: Robert Harkness, December 2007
/
************************************************************************/
 
#include <stdio.h>
#include <mpi.h>
 
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
 
/* function prototypes */
 
int DepositParticleMassField(HierarchyEntry *Grid, FLOAT Time = -1.0);

int CommunicationBufferPurge(void);
 
#ifdef SIB3
int PrepareGravitatingMassField(HierarchyEntry *Grid, int grid1,
				SiblingGridList SiblingList[],
				TopGridData *MetaData, int level,
                                FLOAT When);
#else
int PrepareGravitatingMassField(HierarchyEntry *Grid, TopGridData *MetaData,
                                LevelHierarchyEntry *LevelArray[], int level,
                                FLOAT When);
#endif
 
#ifdef SIB5
int ComputePotentialFieldLevelZero(TopGridData *MetaData,
				   SiblingGridList SiblingList[],
				   HierarchyEntry *Grids[], int NumberOfGrids);
#else
int ComputePotentialFieldLevelZero(TopGridData *MetaData,
				   HierarchyEntry *Grids[], int NumberOfGrids);
#endif

int GenerateGridArray(LevelHierarchyEntry *LevelArray[], int level,
		      HierarchyEntry **Grids[]);
 
 
 
extern int CopyPotentialFieldAverage;
 
#define GRIDS_PER_LOOP 20000
 
// =======================================================================
// This routine sets all the boundary conditions for Grids by either
// interpolating from their parents or copying from sibling grids.
// =======================================================================
 
#ifdef SIB2
int SetBoundaryConditions(HierarchyEntry *Grids[], int NumberOfGrids,
			  SiblingGridList SiblingList[],
			  int level, TopGridData *MetaData,
			  ExternalBoundary *Exterior, LevelHierarchyEntry *Level)
#else
int SetBoundaryConditions(HierarchyEntry *Grids[], int NumberOfGrids,
                          int level, TopGridData *MetaData,
                          ExternalBoundary *Exterior, LevelHierarchyEntry *Level)
#endif
{

  int grid1, grid2;
 
  /* -------------- FIRST PASS ----------------- */

#ifdef FORCE_MSG_PROGRESS 
  MPI_Barrier(MPI_COMM_WORLD);
#endif

  CommunicationDirection = COMMUNICATION_SEND;
 
  if (traceMPI) fprintf(tracePtr, "SBC send\n");

  int *status1 = new int[NumberOfGrids];
  int *status2 = new int[NumberOfGrids];

#pragma omp parallel private(grid1, grid2) shared(level, NumberOfGrids, Grids, MetaData, Exterior, SiblingList, status1, status2) default(none) if(level > 0)
  {
#pragma omp for schedule(dynamic) 
  for (grid1 = 0; grid1 < NumberOfGrids; grid1++) {
 
    /* a) Interpolate boundaries from the parent grid or set external
       boundary conditions. */
 
    if (level == 0) {
      status1[grid1] = Grids[grid1]->GridData->SetExternalBoundaryValues(Exterior);
	// fprintf(stderr, "Error in grid->SetExternalBoundaryValues.\n");
    }
    else {
      status1[grid1] = Grids[grid1]->GridData->InterpolateBoundaryFromParent
	   (Grids[grid1]->ParentGrid->GridData);
	// fprintf(stderr, "Error in grid->InterpolateBoundaryFromParent.\n");
    }

  }

  } // End of parallel region

  for (grid1 = 0; grid1 < NumberOfGrids; grid1++) { 

    /* b) Copy any overlapping zones for sibling grids.  */

#ifdef SIB2
    for (grid2 = 0; grid2 < SiblingList[grid1].NumberOfSiblings; grid2++)
      status2[grid2] = Grids[grid1]->GridData->CheckForOverlap(
				     SiblingList[grid1].GridList[grid2],
				     MetaData->LeftFaceBoundaryCondition,
				     MetaData->RightFaceBoundaryCondition,
				     &grid::CopyZonesFromGrid);
      // fprintf(stderr, "Error in grid->CopyZonesFromGrid.\n");
#else
    for (grid2 = 0; grid2 < NumberOfGrids; grid2++)
      status2[grid2] = Grids[grid1]->GridData->CheckForOverlap(Grids[grid2]->GridData,
                                     MetaData->LeftFaceBoundaryCondition,
                                     MetaData->RightFaceBoundaryCondition,
                                     &grid::CopyZonesFromGrid);
      // fprintf(stderr, "Error in grid->CopyZonesFromGrid.\n");
#endif

    /* c) Apply external reflecting boundary conditions, if needed.  */

    status1[grid1] = Grids[grid1]->GridData->CheckForExternalReflections(
							   MetaData->LeftFaceBoundaryCondition,
							   MetaData->RightFaceBoundaryCondition);
    // fprintf(stderr, "Error in grid->CheckForExternalReflections.\n");
 
  } // end loop over grids


#ifdef FORCE_MSG_PROGRESS 
  MPI_Barrier(MPI_COMM_WORLD);
#endif
 
  /* -------------- SECOND PASS ----------------- */
 
  CommunicationDirection = COMMUNICATION_RECEIVE;
 
  if (traceMPI) fprintf(tracePtr, "SBC recv\n");


#pragma omp parallel private(grid1, grid2) shared(level, NumberOfGrids, Grids, MetaData, SiblingList, status1, status2) default(none) if(level > 0)
  {
#pragma omp for schedule(dynamic)
 
  for (grid1 = 0; grid1 < NumberOfGrids; grid1++) {
 
    /* a) Interpolate boundaries from the parent grid or set external
       boundary conditions. */
 
    if (level > 0)
      status1[grid1] = Grids[grid1]->GridData->InterpolateBoundaryFromParent
	   (Grids[grid1]->ParentGrid->GridData);
	// fprintf(stderr, "Error in grid->InterpolateBoundaryFromParent.\n");

  }

  }  // End of parallel region

  for (grid1 = 0; grid1 < NumberOfGrids; grid1++) {
 
    /* b) Copy any overlapping zones for sibling grids.  */

#ifdef SIB2
    for (grid2 = 0; grid2 < SiblingList[grid1].NumberOfSiblings; grid2++)
      status2[grid2] = Grids[grid1]->GridData->CheckForOverlap(
				     SiblingList[grid1].GridList[grid2],
				     MetaData->LeftFaceBoundaryCondition,
				     MetaData->RightFaceBoundaryCondition,
				     &grid::CopyZonesFromGrid);
      // fprintf(stderr, "Error in grid->CopyZonesFromGrid.\n");
#else
    for (grid2 = 0; grid2 < NumberOfGrids; grid2++)
      status2[grid2] = Grids[grid1]->GridData->CheckForOverlap(Grids[grid2]->GridData,
                                     MetaData->LeftFaceBoundaryCondition,
                                     MetaData->RightFaceBoundaryCondition,
                                     &grid::CopyZonesFromGrid);
      // fprintf(stderr, "Error in grid->CopyZonesFromGrid.\n");
#endif
 
  } // end loop over grids


#ifdef FORCE_MSG_PROGRESS 
  MPI_Barrier(MPI_COMM_WORLD);
#endif

  CommunicationDirection = COMMUNICATION_SEND_RECEIVE;

  delete [] status1;
  delete [] status2;
 
  return SUCCESS;
}
 
 
/* ======================================================================= */
/* This routine prepares the density field for all the grids on this level,
   both particle and baryonic densities.  It also calculates the potential
   field if this is level 0 (since this involves communication). */
 
#ifdef SIB3
int PrepareDensityField(LevelHierarchyEntry *LevelArray[],
			SiblingGridList SiblingList[],
			int level, TopGridData *MetaData, FLOAT When)
#else   // !SIB3
int PrepareDensityField(LevelHierarchyEntry *LevelArray[],
                        int level, TopGridData *MetaData, FLOAT When)
#endif  // end SIB3
{
 
  int grid1, grid2;
 
  /* Set the time for evaluation of the fields, etc. */
 
  FLOAT EvaluateTime = LevelArray[level]->GridData->ReturnTime() +
                   When*LevelArray[level]->GridData->ReturnTimeStep();
 
  /* If level is above MaximumGravityRefinementLevel, then just update the
     gravity at the MaximumGravityRefinementLevel. */
 
  int reallevel = level;
  level = min(level, MaximumGravityRefinementLevel);
 
  /* Create an array (Grids) of all the grids. */
 
  typedef HierarchyEntry* HierarchyEntryPointer;
  HierarchyEntry **Grids;
  int NumberOfGrids = GenerateGridArray(LevelArray, level, &Grids);

  int *status1 = new int[NumberOfGrids];
  int *status2 = new int[NumberOfGrids];
 
  /* Grids: Deposit particles in their GravitatingMassFieldParticles. */
 
  if (traceMPI) fprintf(tracePtr, "PrepareDensityField: Enter DepositParticleMassField (Send)\n");

#ifdef FORCE_MSG_PROGRESS 
  MPI_Barrier(MPI_COMM_WORLD);
#endif

  CommunicationDirection = COMMUNICATION_SEND;

  for (grid1 = 0; grid1 < NumberOfGrids; grid1++)
    if (DepositParticleMassField(Grids[grid1], EvaluateTime) == FAIL) {
      fprintf(stderr, "Error in DepositParticleMassField.\n");
      return FAIL;
    }

#ifdef FORCE_MSG_PROGRESS 
  MPI_Barrier(MPI_COMM_WORLD);
#endif

  if (traceMPI) fprintf(tracePtr, "PrepareDensityField: Enter DepositParticleMassField (Receive)\n");
 
  CommunicationDirection = COMMUNICATION_RECEIVE;

  for (grid1 = 0; grid1 < NumberOfGrids; grid1++)
    if (DepositParticleMassField(Grids[grid1], EvaluateTime) == FAIL) {
      fprintf(stderr, "Error in DepositParticleMassField.\n");
      return FAIL;
    }

#ifdef FORCE_BUFFER_PURGE
  CommunicationBufferPurge();
#endif

#ifdef FORCE_MSG_PROGRESS 
  MPI_Barrier(MPI_COMM_WORLD);
#endif
 
  /* Grids: compute the GravitatingMassField (baryons & particles). */
 
  if (traceMPI) fprintf(tracePtr, "PrepareDensityField: P(%"ISYM"): PGMF1 (send)\n", MyProcessorNumber);
 
  CommunicationDirection = COMMUNICATION_SEND;
 
#ifdef SIB3
  for (grid1 = 0; grid1 < NumberOfGrids; grid1++)
    if (PrepareGravitatingMassField(Grids[grid1], grid1, SiblingList,
				    MetaData, level, When) == FAIL) {
      fprintf(stderr, "Error in PrepareGravitatingMassField.\n");
      return FAIL;
    }
#else
  for (grid1 = 0; grid1 < NumberOfGrids; grid1++)
    if (PrepareGravitatingMassField(Grids[grid1], MetaData, LevelArray,
                                    level, When) == FAIL) {
      fprintf(stderr, "Error in PrepareGravitatingMassField.\n");
      return FAIL;
    }
#endif

#ifdef FORCE_MSG_PROGRESS 
  MPI_Barrier(MPI_COMM_WORLD);
#endif

  if (traceMPI) fprintf(tracePtr, "PrepareDensityField: P(%"ISYM"): PGMF2 (receive)\n", MyProcessorNumber);
 
  CommunicationDirection = COMMUNICATION_RECEIVE;
 
#ifdef SIB3
  for (grid1 = 0; grid1 < NumberOfGrids; grid1++)
    if (PrepareGravitatingMassField(Grids[grid1], grid1, SiblingList,
				    MetaData, level, When) == FAIL) {
      fprintf(stderr, "Error in PrepareGravitatingMassField.\n");
      return FAIL;
    }
#else
  for (grid1 = 0; grid1 < NumberOfGrids; grid1++)
    if (PrepareGravitatingMassField(Grids[grid1], MetaData, LevelArray,
                                    level, When) == FAIL) {
      fprintf(stderr, "Error in PrepareGravitatingMassField.\n");
      return FAIL;
    }
#endif

#ifdef FORCE_BUFFER_PURGE
  CommunicationBufferPurge();
#endif

#ifdef FORCE_MSG_PROGRESS 
  MPI_Barrier(MPI_COMM_WORLD);
#endif
 
  /* Copy overlapping mass fields to ensure consistency and B.C.'s. */
 
  //  if (level > 0)
 
  if (traceMPI) fprintf(tracePtr, "PrepareDensityField: P(%"ISYM"): COMF1 (send)\n", MyProcessorNumber);
 
  CommunicationDirection = COMMUNICATION_SEND;
 
#ifdef SIB1
  for (grid1 = 0; grid1 < NumberOfGrids; grid1++)
    for (grid2 = 0; grid2 < SiblingList[grid1].NumberOfSiblings; grid2++)
      if (Grids[grid1]->GridData->CheckForOverlap(
				   SiblingList[grid1].GridList[grid2],
				   MetaData->LeftFaceBoundaryCondition,
				   MetaData->RightFaceBoundaryCondition,
				   &grid::CopyOverlappingMassField) == FAIL) {
	fprintf(stderr, "Error in grid->CopyOverlappingMassField.\n");
	return FAIL;
      }
#else
  for (grid1 = 0; grid1 < NumberOfGrids; grid1++)
    for (grid2 = 0; grid2 < NumberOfGrids; grid2++)
      if (Grids[grid1]->GridData->CheckForOverlap(Grids[grid2]->GridData,
                                   MetaData->LeftFaceBoundaryCondition,
                                   MetaData->RightFaceBoundaryCondition,
                                   &grid::CopyOverlappingMassField) == FAIL) {
        fprintf(stderr, "Error in grid->CopyOverlappingMassField.\n");
        return FAIL;
      }
#endif

#ifdef FORCE_MSG_PROGRESS 
  MPI_Barrier(MPI_COMM_WORLD);
#endif

  if (traceMPI) fprintf(tracePtr, "PrepareDensityField: P(%"ISYM"): COMF2 (receive)\n", MyProcessorNumber);
 
  CommunicationDirection = COMMUNICATION_RECEIVE;
 
#ifdef SIB1
  for (grid1 = 0; grid1 < NumberOfGrids; grid1++)
    for (grid2 = 0; grid2 < SiblingList[grid1].NumberOfSiblings; grid2++)
      if (Grids[grid1]->GridData->CheckForOverlap(
				   SiblingList[grid1].GridList[grid2],
				   MetaData->LeftFaceBoundaryCondition,
				   MetaData->RightFaceBoundaryCondition,
				   &grid::CopyOverlappingMassField) == FAIL) {
	fprintf(stderr, "Error in grid->CopyOverlappingMassField.\n");
	return FAIL;
      }
#else
  for (grid1 = 0; grid1 < NumberOfGrids; grid1++)
    for (grid2 = 0; grid2 < NumberOfGrids; grid2++)
      if (Grids[grid1]->GridData->CheckForOverlap(Grids[grid2]->GridData,
                                   MetaData->LeftFaceBoundaryCondition,
                                   MetaData->RightFaceBoundaryCondition,
                                   &grid::CopyOverlappingMassField) == FAIL) {
        fprintf(stderr, "Error in grid->CopyOverlappingMassField.\n");
        return FAIL;
      }
#endif
 
#ifdef FORCE_BUFFER_PURGE
  CommunicationBufferPurge();
#endif

#ifdef FORCE_MSG_PROGRESS 
  MPI_Barrier(MPI_COMM_WORLD);
#endif

  CommunicationDirection = COMMUNICATION_SEND_RECEIVE;
 
  /* Compute the potential for the top grid. */
 
  if (level == 0) {
    if (traceMPI) fprintf(tracePtr, "PrepareDensityField: P(%"ISYM"): CPFLZero (send-receive)\n", MyProcessorNumber);

#ifdef SIB5
    if (ComputePotentialFieldLevelZero(MetaData, SiblingList, Grids, NumberOfGrids) == FAIL) {
#else
    if (ComputePotentialFieldLevelZero(MetaData, Grids, NumberOfGrids) == FAIL) {
#endif
      fprintf(stderr, "Error in ComputePotentialFieldLevelZero.\n");
      return FAIL;
    }
  }
 
  /* Compute a first iteration of the potential and share BV's. */
 
#define ITERATE_POTENTIAL
#ifdef ITERATE_POTENTIAL
      if (level > 0) {
	CopyPotentialFieldAverage = 1;
	for (int iterate = 0; iterate < PotentialIterations; iterate++) {
 
	  if (iterate > 0)
	    CopyPotentialFieldAverage = 2;
 
	  int Dummy;
#pragma omp parallel private(grid1, Dummy) shared(NumberOfGrids, Grids, level, CopyGravPotential, EvaluateTime, status1, status2) default(none)
  {
#pragma omp for schedule(dynamic)
	  for (grid1 = 0; grid1 < NumberOfGrids; grid1++) {

	    status1[grid1] = Grids[grid1]->GridData->SolveForPotential(Dummy, level, EvaluateTime);

            if (CopyGravPotential) {
              status2[grid1] = Grids[grid1]->GridData->CopyPotentialToBaryonField();
            }

          } // End of loop over grids

  } // End parallel region
 
          if (traceMPI) fprintf(tracePtr, "ITPOT send\n");

#ifdef FORCE_MSG_PROGRESS 
          MPI_Barrier(MPI_COMM_WORLD);
#endif

	  CommunicationDirection = COMMUNICATION_SEND;
 
#ifdef SIB4
	  for (grid1 = 0; grid1 < NumberOfGrids; grid1++) {
 
            // fprintf(stderr, "#SIBSend on cpu %"ISYM": %"ISYM"\n", MyProcessorNumber, SiblingList[grid1].NumberOfSiblings);
 
            // for (grid2 = SiblingList[grid1].NumberOfSiblings-1; grid2 = 0; grid2--)
	    for (grid2 = 0; grid2 < SiblingList[grid1].NumberOfSiblings; grid2++)
	     if (Grids[grid1]->GridData->CheckForOverlap(
				   SiblingList[grid1].GridList[grid2],
				   MetaData->LeftFaceBoundaryCondition,
				   MetaData->RightFaceBoundaryCondition,
				   &grid::CopyPotentialField) == FAIL) {
	       fprintf(stderr, "Error in grid->CopyPotentialField.\n");
	       return FAIL;
	     }
 
            grid2 = grid1;
             if (Grids[grid1]->GridData->CheckForOverlap(
                                   Grids[grid2]->GridData,
                                   MetaData->LeftFaceBoundaryCondition,
                                   MetaData->RightFaceBoundaryCondition,
                                   &grid::CopyPotentialField) == FAIL) {
               fprintf(stderr, "Error in grid->CopyPotentialField.\n");
               return FAIL;
             }
 
           }
#else
          for (grid1 = 0; grid1 < NumberOfGrids; grid1++)
            for (grid2 = 0; grid2 < NumberOfGrids; grid2++)
             if (Grids[grid1]->GridData->CheckForOverlap(Grids[grid2]->GridData,
                                   MetaData->LeftFaceBoundaryCondition,
                                   MetaData->RightFaceBoundaryCondition,
                                   &grid::CopyPotentialField) == FAIL) {
               fprintf(stderr, "Error in grid->CopyPotentialField.\n");
               return FAIL;
             }
#endif

#ifdef FORCE_MSG_PROGRESS 
          MPI_Barrier(MPI_COMM_WORLD);
#endif

          if (traceMPI) fprintf(tracePtr, "ITPOT recv\n");
 
	  CommunicationDirection = COMMUNICATION_RECEIVE;
 
#ifdef SIB4
	  for (grid1 = 0; grid1 < NumberOfGrids; grid1++) {
 
            // fprintf(stderr, "#SIBRecv on cpu %"ISYM": %"ISYM"\n", MyProcessorNumber, SiblingList[grid1].NumberOfSiblings);
 
            // for (grid2 = SiblingList[grid1].NumberOfSiblings-1; grid2 = 0; grid2--)
	    for (grid2 = 0; grid2 < SiblingList[grid1].NumberOfSiblings; grid2++)
	     if (Grids[grid1]->GridData->CheckForOverlap(
				   SiblingList[grid1].GridList[grid2],
				   MetaData->LeftFaceBoundaryCondition,
				   MetaData->RightFaceBoundaryCondition,
				   &grid::CopyPotentialField) == FAIL) {
	       fprintf(stderr, "Error in grid->CopyPotentialField.\n");
	       return FAIL;
	     }
 
            grid2 = grid1;
             if (Grids[grid1]->GridData->CheckForOverlap(
                                   Grids[grid2]->GridData,
                                   MetaData->LeftFaceBoundaryCondition,
                                   MetaData->RightFaceBoundaryCondition,
                                   &grid::CopyPotentialField) == FAIL) {
               fprintf(stderr, "Error in grid->CopyPotentialField.\n");
               return FAIL;
             }
 
           }
#else
          for (grid1 = 0; grid1 < NumberOfGrids; grid1++)
            for (grid2 = 0; grid2 < NumberOfGrids; grid2++)
             if (Grids[grid1]->GridData->CheckForOverlap(Grids[grid2]->GridData,
                                   MetaData->LeftFaceBoundaryCondition,
                                   MetaData->RightFaceBoundaryCondition,
                                   &grid::CopyPotentialField) == FAIL) {
               fprintf(stderr, "Error in grid->CopyPotentialField.\n");
               return FAIL;
             }
#endif

#ifdef FORCE_MSG_PROGRESS 
          MPI_Barrier(MPI_COMM_WORLD);
#endif

	  CommunicationDirection = COMMUNICATION_SEND_RECEIVE;
 
 
	}
	CopyPotentialFieldAverage = 0;
      }
#endif /* ITERATE_POTENTIAL */
 
  /* if level > MaximumGravityRefinementLevel, then do final potential
     solve (and acceleration interpolation) here rather than in the main
     EvolveLevel since it involves communications. */
 
  if (reallevel > MaximumGravityRefinementLevel) {
 
    /* compute potential and acceleration on coarser level [LOCAL]
       (but only if there is at least a subgrid -- it should be only
        if there is a subgrrid on reallevel, but this is ok). */
 
    for (grid1 = 0; grid1 < NumberOfGrids; grid1++)
      if (Grids[grid1]->NextGridNextLevel != NULL) {
	Grids[grid1]->GridData->SolveForPotential(level,
					       MaximumGravityRefinementLevel);
        if (CopyGravPotential)
          Grids[grid1]->GridData->CopyPotentialToBaryonField();
        else
  	  Grids[grid1]->GridData->ComputeAccelerationField(
                           (HydroMethod == Zeus_Hydro) ? DIFFERENCE_TYPE_STAGGERED : DIFFERENCE_TYPE_NORMAL,
					       MaximumGravityRefinementLevel);
      }
 
    /* Interpolate potential for reallevel grids from coarser grids. */
 
    if (!CopyGravPotential) {
 
      int Dummy;
      LevelHierarchyEntry *Temp = LevelArray[reallevel];

#ifdef FORCE_MSG_PROGRESS 
      MPI_Barrier(MPI_COMM_WORLD);
#endif

      CommunicationDirection = COMMUNICATION_SEND;

      while (Temp != NULL) {
        HierarchyEntry *Temp3 = Temp->GridHierarchyEntry;
        for (Dummy = reallevel; Dummy > MaximumGravityRefinementLevel; Dummy--)
  	  Temp3 = Temp3->ParentGrid;
        if (Temp->GridData->InterpolateAccelerations(Temp3->GridData) == FAIL) {
	  fprintf(stderr, "Error in grid->InterpolateAccelerations.\n");
	  return FAIL;
        }
        Temp = Temp->NextGridThisLevel;
      }

#ifdef FORCE_MSG_PROGRESS 
      MPI_Barrier(MPI_COMM_WORLD);
#endif

      CommunicationDirection = COMMUNICATION_RECEIVE;

      Temp = LevelArray[reallevel];

      while (Temp != NULL) {
        HierarchyEntry *Temp3 = Temp->GridHierarchyEntry;
        for (Dummy = reallevel; Dummy > MaximumGravityRefinementLevel; Dummy--)
	  Temp3 = Temp3->ParentGrid;
        if (Temp->GridData->InterpolateAccelerations(Temp3->GridData) == FAIL) {
	  fprintf(stderr, "Error in grid->InterpolateAccelerations.\n");
	  return FAIL;
        }
        Temp = Temp->NextGridThisLevel;
      }

#ifdef FORCE_MSG_PROGRESS 
      MPI_Barrier(MPI_COMM_WORLD);
#endif

      CommunicationDirection = COMMUNICATION_SEND_RECEIVE;

    } // end:  if (!CopyGravPotential)
 
  } // end: if (reallevel > MaximumGravityRefinementLevel)

  // --------------------------------------------------
  // MEMORY LEAK FIX
  //
  // Adding missing delete [] () for Grids[] allocated in
  // GenerateGridArray()
  // --------------------------------------------------

  delete [] Grids;

  // --------------------------------------------------

  delete [] status1;
  delete [] status2;
 
  return SUCCESS;
}
 
 
/* ======================================================================= */
/* This routines does the flux correction and project for all grids on this
   level from the list of subgrids. */
 
#ifdef FLUX_FIX
int UpdateFromFinerGrids(int level, HierarchyEntry *Grids[], int NumberOfGrids,
			 int NumberOfSubgrids[],
			 fluxes **SubgridFluxesEstimate[],
			 LevelHierarchyEntry* SUBlingList[],
			 TopGridData *MetaData)
#else
int UpdateFromFinerGrids(int level, HierarchyEntry *Grids[], int NumberOfGrids,
			 int NumberOfSubgrids[],
			 fluxes **SubgridFluxesEstimate[])
#endif
 
{
 
  int grid1, subgrid;
  HierarchyEntry *NextGrid;
 
#ifdef FLUX_FIX
  int SUBlingGrid;
  LevelHierarchyEntry *NextEntry;
#endif
 
  /* Define a temporary flux holder for the refined fluxes. */
 
  fluxes SubgridFluxesRefined;
 
  /* For each grid,
     (a) project the subgrid's solution into this grid (step #18)
     (b) correct for the difference between this grid's fluxes and the
         subgrid's fluxes. (step #19) */
 
  /* -------------- FIRST PASS ----------------- */

#ifdef FORCE_MSG_PROGRESS 
  MPI_Barrier(MPI_COMM_WORLD);
#endif

  CommunicationDirection = COMMUNICATION_SEND;
 
  for (grid1 = 0; grid1 < NumberOfGrids; grid1++) {
 
    /* Loop over subgrids for this grid. */
 
    NextGrid = Grids[grid1]->NextGridNextLevel;
    subgrid = 0;
    while (NextGrid != NULL && FluxCorrection) {
 
      /* Project subgrid's refined fluxes to the level of this grid. */
 
      if (NextGrid->GridData->GetProjectedBoundaryFluxes(
		      Grids[grid1]->GridData, SubgridFluxesRefined) == FAIL) {
	fprintf(stderr, "Error in grid->GetProjectedBoundaryFluxes.\n");
	return FAIL;
      }
 
      NextGrid = NextGrid->NextGridThisLevel;
      subgrid++;
    }
 
#ifdef FLUX_FIX
    NextEntry = SUBlingList[grid1];
 
    while (NextEntry != NULL && FluxCorrection) {
      /* make sure this isn't a "proper" subgrid */
      if( !(NextEntry->GridHierarchyEntry->ParentGrid == Grids[grid1]) ){
        /* Project subgrid's refined fluxes to the level of this grid. */
        if (NextEntry->GridData->
            GetProjectedBoundaryFluxes( Grids[grid1]->GridData,
                                        SubgridFluxesRefined ) == FAIL) {
          fprintf(stderr, "Error in grid->GetProjectedBoundaryFluxes.\n");
          return FAIL;
        }
      }
      NextEntry = NextEntry->NextGridThisLevel;
    }
#endif
 
    /* Loop over subgrids for this grid: replace solution. */
 
    NextGrid = Grids[grid1]->NextGridNextLevel;
    while (NextGrid != NULL) {
 
      /* Project the subgrid solution into this grid. */
 
      if (NextGrid->GridData->ProjectSolutionToParentGrid
	                                   (*Grids[grid1]->GridData) == FAIL) {
	fprintf(stderr, "Error in grid->ProjectSolutionToParentGrid.\n");
	return FAIL;
      }
 
      NextGrid = NextGrid->NextGridThisLevel;
    }
 
  } // end of loop over subgrids

#ifdef FORCE_MSG_PROGRESS 
  MPI_Barrier(MPI_COMM_WORLD);
#endif
 
  /* -------------- SECOND PASS ----------------- */
 
  CommunicationDirection = COMMUNICATION_RECEIVE;
 
  for (grid1 = 0; grid1 < NumberOfGrids; grid1++) {
 
    /* Loop over subgrids for this grid. */
    /* May have to do it this way to ensure that the fluxes are
       matched up to the correct grids. */
 
    NextGrid = Grids[grid1]->NextGridNextLevel;
    subgrid = 0;
    while (NextGrid != NULL && FluxCorrection) {
 
      /* Project subgrid's refined fluxes to the level of this grid. */
 
      if (NextGrid->GridData->GetProjectedBoundaryFluxes(
		      Grids[grid1]->GridData, SubgridFluxesRefined) == FAIL) {
	fprintf(stderr, "Error in grid->GetProjectedBoundaryFluxes.\n");
	return FAIL;
      }
	
      /* Correct this grid for the refined fluxes (step #19)
	 (this also deletes the fields in SubgridFluxesRefined). */
 
#ifdef FLUX_FIX
      if (Grids[grid1]->GridData->CorrectForRefinedFluxes
          (SubgridFluxesEstimate[grid1][subgrid], &SubgridFluxesRefined,
           SubgridFluxesEstimate[grid1][NumberOfSubgrids[grid1] - 1],
           FALSE, MetaData)
          == FAIL) {
        fprintf(stderr, "Error in grid->CorrectForRefinedFluxes.\n");
        return FAIL;
      }
#else
      if (Grids[grid1]->GridData->CorrectForRefinedFluxes
	  (SubgridFluxesEstimate[grid1][subgrid], &SubgridFluxesRefined,
	   SubgridFluxesEstimate[grid1][NumberOfSubgrids[grid1] - 1]     )
	  == FAIL) {
	fprintf(stderr, "Error in grid->CorrectForRefinedFluxes.\n");
	return FAIL;
      }
#endif
 
      NextGrid = NextGrid->NextGridThisLevel;
      subgrid++;
    }
 
#ifdef FLUX_FIX /* Do flux corrections from Sublings */
 
    NextEntry = SUBlingList[grid1];
 
    while (NextEntry != NULL && FluxCorrection) {
      /* make sure this isn't a "proper" subgrid */
      if( NextEntry->GridHierarchyEntry->ParentGrid != Grids[grid1] ){
 
        /* Project subgrid's refined fluxes to the level of this grid. */
        if (NextEntry->GridData->
            GetProjectedBoundaryFluxes( Grids[grid1]->GridData,
                                        SubgridFluxesRefined ) == FAIL) {
          fprintf(stderr, "Error in grid->GetProjectedBoundaryFluxes.\n");
          return FAIL;
        }
 
        /* Correct this grid for the refined fluxes (step #19)
           (this also deletes the fields in SubgridFluxesRefined). */
 
        if (Grids[grid1]->GridData->CorrectForRefinedFluxes
            (SubgridFluxesEstimate[grid1][NumberOfSubgrids[grid1] - 1],
             &SubgridFluxesRefined,
             SubgridFluxesEstimate[grid1][NumberOfSubgrids[grid1] - 1],
             TRUE, MetaData) == FAIL) {
          fprintf(stderr, "Error in grid->CorrectForRefinedFluxes.\n");
          return FAIL;
        }
      }
      NextEntry = NextEntry->NextGridThisLevel;
    }
#endif
 
    /* Loop over subgrids for this grid: replace solution. */
 
    NextGrid = Grids[grid1]->NextGridNextLevel;
 
    while (NextGrid != NULL) {
 
      /* Project the subgrid solution into this grid. */
 
      if (NextGrid->GridData->ProjectSolutionToParentGrid
	                                   (*Grids[grid1]->GridData) == FAIL) {
	fprintf(stderr, "Error in grid->ProjectSolutionToParentGrid.\n");
	return FAIL;
      }
 
      NextGrid = NextGrid->NextGridThisLevel;
    }
 
  } // end of loop over subgrids

#ifdef FORCE_MSG_PROGRESS 
  MPI_Barrier(MPI_COMM_WORLD);
#endif

  CommunicationDirection = COMMUNICATION_SEND_RECEIVE;
 
  return SUCCESS;
}
 
 
 
/* ======================================================================= */
/* This routine simply converts a linked list of grids into an array of
   pointers. */
 
int GenerateGridArray(LevelHierarchyEntry *LevelArray[], int level,
		      HierarchyEntry **Grids[])
{
 
  /* Count the number of grids on this level. */
 
  int NumberOfGrids = 0, counter = 0;
  LevelHierarchyEntry *Temp = LevelArray[level];
  while (Temp != NULL) {
    NumberOfGrids++;
    Temp             = Temp->NextGridThisLevel;
  }
 
  /* Create a list of pointers and number of subgrids (and fill it out). */
 
  typedef HierarchyEntry* HierarchyEntryPointer;
  *Grids = new HierarchyEntryPointer[NumberOfGrids];
  Temp = LevelArray[level];
  while (Temp != NULL) {
    (*Grids)[counter++] = Temp->GridHierarchyEntry;
    Temp              = Temp->NextGridThisLevel;
  }
 
  return NumberOfGrids;
}
 
