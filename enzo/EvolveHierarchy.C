/***********************************************************************
/
/  EVOLVE HIERARCHY FUNCTION
/
/  written by: Greg Bryan
/  date:       November, 1994
/  modified1:  February, 1995 by GB
/              Changed to reflect changes in EvolveGrid & EvolveTopGrid.
/  modified2:  July, 1995 by GB
/              Changed to reflect new routine EvolveLevel.
/  modified3:  February, 2006 by Daniel Reynolds
/              Updated call interface to ComputePotentialFieldLevelZero
/  modified4:  February, 2007 by Robert Harkness
/              Group and in-core i/o
/  modified5:  December, 2007 by Robert Harkness
/              Remove calls to ComputePotential to allow use of the
/              Fast Sibling Locator to remedy Ncpu^2 scaling
/  modified6:  February, 2008 by Robert Harkness
/              Conditional calls to MPI_Barrier to force MPICH progress
/
/  PURPOSE:
/    This routine is responsible for the evolution of the grid hierarchy.
/    It assumes the hierarchy is already constructed and the grids
/    initialized.  The routine then loops over time until one of the
/    stopping criteria is reached.  This routine also handles data dumps,
/    history dumps and restarts dumps (although the later two have not
/    yet been implemented).
/
************************************************************************/

#ifdef RAD_HYDRO
#include "ImplicitProblemABC_preincludes.h"
#endif
 
#include <stdio.h>
 
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
#include "CosmologyParameters.h"
#ifdef RAD_HYDRO
#include "ImplicitProblemABC.h"
#endif
 
// function prototypes
 
int RebuildHierarchy(TopGridData *MetaData,
		     LevelHierarchyEntry *LevelArray[], int level);

int EvolveLevel(TopGridData *MetaData, LevelHierarchyEntry *LevelArray[],
		int level, float dtLevelAbove, ExternalBoundary *Exterior
#ifdef RAD_HYDRO
		, ImplicitProblemABC *ImplicitSolver
#endif
);

int WriteAllData(char *basename, int filenumber,
		 HierarchyEntry *TopGrid, TopGridData &MetaData,
		 ExternalBoundary *Exterior,
#ifdef RAD_HYDRO
		 ImplicitProblemABC *ImplicitSolver,
#endif
		 FLOAT WriteTime = -1);

int Group_WriteAllData(char *basename, int filenumber,
		 HierarchyEntry *TopGrid, TopGridData &MetaData,
		 ExternalBoundary *Exterior,
#ifdef RAD_HYDRO
		 ImplicitProblemABC *ImplicitSolver,
#endif
		 FLOAT WriteTime = -1);

int CopyOverlappingZones(grid* CurrentGrid, TopGridData *MetaData,
			 LevelHierarchyEntry *LevelArray[], int level);
int TestGravityCheckResults(LevelHierarchyEntry *LevelArray[]);
int TestGravitySphereCheckResults(LevelHierarchyEntry *LevelArray[]);
int CheckForOutput(HierarchyEntry *TopGrid, TopGridData &MetaData,
		   ExternalBoundary *Exterior,
#ifdef RAD_HYDRO
		   ImplicitProblemABC *ImplicitSolver,
#endif
		   int &WroteData);
int CheckForTimeAction(LevelHierarchyEntry *LevelArray[],
		       TopGridData &MetaData);
int CosmologyComputeExpansionFactor(FLOAT time, FLOAT *a, FLOAT *dadt);
// int CountGrids(TopGridData &MetaData, LevelHierarchyEntry *LevelArray[]);
int OutputLevelInformation(FILE *fptr, TopGridData &MetaData,
			   LevelHierarchyEntry *LevelArray[]);
int PrepareGravitatingMassField(HierarchyEntry *Grid, TopGridData *MetaData,
				LevelHierarchyEntry *LevelArray[], int level);
float CommunicationMinValue(float Value);
int ReduceFragmentation(HierarchyEntry &TopGrid, TopGridData &MetaData,
			ExternalBoundary *Exterior,
			LevelHierarchyEntry *LevelArray[]);

int SetDumpFlagsOn(void);
int SetDumpFlagsMask(void);
int CheckForHalt(void);
int CheckForDump(void);
int CheckForCycleDump(void);

#ifdef MEM_TRACE
Eint64 mused(void);
#endif
 
 
#define NO_REDUCE_FRAGMENTATION
 
 
 
 
int EvolveHierarchy(HierarchyEntry &TopGrid, TopGridData &MetaData,
                    ExternalBoundary *Exterior,
#ifdef RAD_HYDRO
		    ImplicitProblemABC *ImplicitSolver,
#endif
		    LevelHierarchyEntry *LevelArray[], float Initialdt)
{
 
  float dt;
 
  int i, Stop = FALSE, WroteData, HaltFlag, DumpFlag, DumpStep;
  double tlev0, tlev1, treb0, treb1, tloop0, tloop1, tentry, texit;

#ifdef USE_MPI
  tentry = MPI_Wtime();
#endif
 
  if (MetaData.Time        >= MetaData.StopTime ) Stop = TRUE;
  if (MetaData.CycleNumber >= MetaData.StopCycle) Stop = TRUE;
 
#ifdef USE_MPI
  MetaData.CPUTime = MPI_Wtime();
#endif

#ifdef MEM_TRACE
    Eint64 MemInUse;
#endif

  SetDumpFlagsOn(); 

  /* Attach RandomForcingFields to BaryonFields temporarily to apply BCs */
 
  if (RandomForcing) { //AK
    LevelHierarchyEntry *Temp = LevelArray[0];
    while (Temp != NULL) {
      Temp->GridData->AppendForcingToBaryonFields();
      Temp = Temp->NextGridThisLevel;
    }
    Exterior->AppendForcingToBaryonFields();
  }
 
  /* Set top grid boundary conditions. */
 
  LevelHierarchyEntry *Temp = LevelArray[0];

#ifdef MEM_TRACE
    MemInUse = mused();
    fprintf(memtracePtr, "Enter EH %8"ISYM"  %16"ISYM" \n", MetaData.CycleNumber, MemInUse);
    fflush(memtracePtr);
#endif

#ifdef FORCE_MSG_PROGRESS
  MPI_Barrier(MPI_COMM_WORLD);
#endif

  CommunicationDirection = COMMUNICATION_SEND;

  while (Temp != NULL) {
    if (Temp->GridData->SetExternalBoundaryValues(Exterior) == FAIL) {
      fprintf(stderr, "Error in grid->SetExternalBoundaryValues.\n");
      return FAIL;
    }
    if (CopyOverlappingZones(Temp->GridData, &MetaData, LevelArray, 0) == FAIL) {
      fprintf(stderr, "Error in CopyOverlappingZones.\n");
      return FAIL;
    }
    Temp = Temp->NextGridThisLevel;
  }

#ifdef FORCE_MSG_PROGRESS 
  MPI_Barrier(MPI_COMM_WORLD);
#endif

  CommunicationDirection = COMMUNICATION_RECEIVE;

  Temp = LevelArray[0];
  while (Temp != NULL) {
    if (Temp->GridData->SetExternalBoundaryValues(Exterior) == FAIL) {
      fprintf(stderr, "Error in grid->SetExternalBoundaryValues.\n");
      return FAIL;
    }
    if (CopyOverlappingZones(Temp->GridData, &MetaData, LevelArray, 0) == FAIL) {
      fprintf(stderr, "Error in CopyOverlappingZones.\n");
      return FAIL;
    }
    Temp = Temp->NextGridThisLevel;
  }
 
  CommunicationDirection = COMMUNICATION_SEND_RECEIVE;

#ifdef FORCE_MSG_PROGRESS
  MPI_Barrier(MPI_COMM_WORLD);
#endif
 
#ifdef MEM_TRACE
    MemInUse = mused();
    fprintf(memtracePtr, "Bdry set %8"ISYM"  %16"ISYM" \n", MetaData.CycleNumber, MemInUse);
    fflush(memtracePtr);
#endif
 
  /* Remove RandomForcingFields from BaryonFields when BCs are set. */
 
  if (RandomForcing) { //AK
    LevelHierarchyEntry *Temp = LevelArray[0];
    while (Temp != NULL) {
      Temp->GridData->DetachForcingFromBaryonFields();
      Temp = Temp->NextGridThisLevel;
    }
    Exterior->DetachForcingFromBaryonFields();
  }
 
  /* Check for output. */
 
  if (CheckForOutput(&TopGrid, MetaData, Exterior,
#ifdef RAD_HYDRO
		     ImplicitSolver,
#endif
		     WroteData) == FAIL) {
    fprintf(stderr, "Error in CheckForOutput.\n");
    return FAIL;
  }

#ifdef MEM_TRACE
    MemInUse = mused();
    fprintf(memtracePtr, "Output %8"ISYM"  %16"ISYM" \n", MetaData.CycleNumber, MemInUse);
    fflush(memtracePtr);
#endif
 
  /* Compute the acceleration field so ComputeTimeStep can find dtAccel.
     (Actually, this is a huge pain-in-the-ass, so only do it if the
      problem really requires it). */
 
/*
  if (ProblemType == 21) {
    PrepareGravitatingMassField(&TopGrid, &MetaData, LevelArray, 0);
    ComputePotentialFieldLevelZero(&MetaData, Grids, NumberOfGrids);
    TopGrid.GridData->ComputeAccelerationField(GRIDS);
  }
*/
 
  /* Do the first grid regeneration. */

  // CountGrids(MetaData, LevelArray);
 
  if (RebuildHierarchy(&MetaData, LevelArray, 0) == FAIL) {
    fprintf(stderr, "Error in RebuildHierarchy.\n");
    return FAIL;
  }

#ifdef MEM_TRACE
    MemInUse = mused();
    fprintf(memtracePtr, "1st rebuild %8"ISYM"  %16"ISYM" \n", MetaData.CycleNumber, MemInUse);
    fflush(memtracePtr);
#endif
 
  /* Open the OutputLevelInformation file. */
 
  FILE *LevelInfofptr;
 
  if (MyProcessorNumber == ROOT_PROCESSOR) {
    LevelInfofptr = fopen("OutputLevelInformation.out", "w");
    fclose(LevelInfofptr);
  }
 
#ifdef USE_JBPERF
  Eint32 jb_iter;
#endif

  /* ====== MAIN LOOP ===== */

  while (!Stop) {

#ifdef USE_JBPERF
    jb_iter = MetaData.CycleNumber;
    static bool isFirstCall = true;
    if ((jb_iter % JB_ITER_PER_SEGMENT)==0 || isFirstCall) jbPerf.begin("EL");
    isFirstCall = false;
    jbPerf.attribute ("timestep",&jb_iter, JB_INT);
    jbPerf.start("EL");
    int time_sim = 1000000*MetaData.Time;
    jbPerf.assign("time-sim",time_sim);
#endif

#ifdef USE_MPI
    tloop0 = MPI_Wtime();
#endif

#ifdef MEM_TRACE
    MemInUse = mused();
    fprintf(memtracePtr, "Top %8"ISYM"  %16"ISYM" \n", MetaData.CycleNumber, MemInUse);
    fflush(memtracePtr);
#endif
 
    /* Output level information to log file. */
 
    if (MyProcessorNumber == ROOT_PROCESSOR) {
      LevelInfofptr = fopen("OutputLevelInformation.out", "a");
    }

    // OutputLevelInformation() only needs to be called by all processors
    // when jbPerf is enabled.

    OutputLevelInformation(LevelInfofptr, MetaData, LevelArray);

    // CountGrids(MetaData, LevelArray);

    if (MyProcessorNumber == ROOT_PROCESSOR) {
      fclose(LevelInfofptr);
    }
 
    /* Compute minimum timestep on the top level. */
 
    float dtProc   = huge_number;
    Temp = LevelArray[0];
 
    while (Temp != NULL) {
      dtProc = min(dtProc, Temp->GridData->ComputeTimeStep());
      Temp = Temp->NextGridThisLevel;
    }

    dt = RootGridCourantSafetyNumber*CommunicationMinValue(dtProc);
 
    if (Initialdt != 0) {
      dt = min(dt, Initialdt);
      Initialdt = 0;
    }
 
    /* Make sure timestep doesn't go past an output. */
 
    if (ComovingCoordinates)
      for (i = 0; i < MAX_NUMBER_OF_OUTPUT_REDSHIFTS; i++)
	if (CosmologyOutputRedshift[i] != -1)
	  dt = min(1.0001*(CosmologyOutputRedshiftTime[i]-MetaData.Time), dt);

    for (i = 0; i < MAX_TIME_ACTIONS; i++)
      if (TimeActionTime[i] > 0 && TimeActionType[i] > 0)
	dt = min(1.0001*(TimeActionTime[i] - MetaData.Time), dt);

    if (MetaData.dtDataDump > 0.0) {

      while (MetaData.TimeLastDataDump+MetaData.dtDataDump < MetaData.Time)
	MetaData.TimeLastDataDump += MetaData.dtDataDump;

      dt = min(1.0001*(MetaData.TimeLastDataDump + MetaData.dtDataDump -
		       MetaData.Time), dt);
    }
 
    /* Set the time step.  If it will cause Time += dt > StopTime, then
       set dt = StopTime - Time */
 
    dt = min(MetaData.StopTime - MetaData.Time, dt);
    Temp = LevelArray[0];
 
    while (Temp != NULL) {
      Temp->GridData->SetTimeStep(dt);
      Temp = Temp->NextGridThisLevel;
    }
 
    if (debug) {
      printf("TopGrid dt = %"FSYM"     time = %"GOUTSYM"    cycle = %"ISYM,
	     dt, MetaData.Time, MetaData.CycleNumber);
      if (ComovingCoordinates) {
	FLOAT a, dadt;
	CosmologyComputeExpansionFactor(MetaData.Time, &a, &dadt);
	printf("    z = %"GOUTSYM, (1 + InitialRedshift)/a - 1);
      }
      printf("\n");
    }
 
    /* Evolve the top grid (and hence the entire hierarchy). */

#ifdef USE_MPI 
    MPI_Barrier(MPI_COMM_WORLD);
    tlev0 = MPI_Wtime();
#endif

    // Zero the level timers

    for(i = 0; i < MAX_DEPTH_OF_HIERARCHY; i++) {
      rebuild_timer[i] = 0.0;
      level_timer[i] = 0.0;
      level_count[i] = 0;
    }

    /* Zeroing out the rootgrid Emissivity before EvolveLevel is called
       so when rootgrid emissivity values are calculated they are put in
       clean rootgrid array */
#ifdef EMISSIVITY
    if(StarMakerEmissivityField > 0){
      LevelHierarchyEntry *RootTemp;
      RootTemp = LevelArray[0];
      while (RootTemp != NULL) {
        RootTemp->GridData->ClearEmissivity();
        RootTemp = RootTemp->NextGridThisLevel;
      }
    }
#endif

#ifdef RAD_HYDRO
    if (EvolveLevel(&MetaData, LevelArray, 0, dt, Exterior, ImplicitSolver) == FAIL) {
      if (MyProcessorNumber == ROOT_PROCESSOR) {
	fprintf(stderr, "Error in EvolveLevel.\n");
      }
      return FAIL;
    }
#else
    if (EvolveLevel(&MetaData, LevelArray, 0, dt, Exterior) == FAIL) {
      if (MyProcessorNumber == ROOT_PROCESSOR) {
        fprintf(stderr, "Error in EvolveLevel.\n");
      }
      return FAIL;
    }
#endif


    // Special call for single field dump at high frequency

#ifdef HF_FIELD_DUMP
    if (MetaData.CycleNumber % 5 == 0 ) {

      SetDumpFlagsMask();

#ifdef USE_HDF5_GROUPS

#ifdef RAD_HYDRO
      Group_WriteAllData(MetaData.HistoryDumpName, MetaData.CycleNumber, &TopGrid, MetaData, Exterior, ImplicitSolver);
#else
      Group_WriteAllData(MetaData.HistoryDumpName, MetaData.CycleNumber, &TopGrid, MetaData, Exterior);
#endif 
      
#else

#ifdef RAD_HYDRO
      WriteAllData(MetaData.HistoryDumpName, MetaData.CycleNumber, &TopGrid, MetaData, Exterior, ImplicitSolver);
#else
      WriteAllData(MetaData.HistoryDumpName, MetaData.CycleNumber, &TopGrid, MetaData, Exterior);
#endif

#endif

      SetDumpFlagsOn();

    }
#endif

    DumpFlag = CheckForDump();
    if (DumpFlag) {

      if (MyProcessorNumber == ROOT_PROCESSOR) printf("Checkpoint on DUMP flag.\n");

#ifdef USE_HDF5_GROUPS

#ifdef RAD_HYDRO
      Group_WriteAllData(MetaData.HistoryDumpName, MetaData.CycleNumber, &TopGrid, MetaData, Exterior, ImplicitSolver);
#else
      Group_WriteAllData(MetaData.HistoryDumpName, MetaData.CycleNumber, &TopGrid, MetaData, Exterior);
#endif 
      
#else

#ifdef RAD_HYDRO
      WriteAllData(MetaData.HistoryDumpName, MetaData.CycleNumber, &TopGrid, MetaData, Exterior, ImplicitSolver);
#else
      WriteAllData(MetaData.HistoryDumpName, MetaData.CycleNumber, &TopGrid, MetaData, Exterior);
#endif

#endif

    }

    // Dump on cycles

#ifdef HF_CYCLE_DUMP
    DumpStep = CheckForCycleDump();

    if (MetaData.CycleNumber % DumpStep == 0) {

      if (MyProcessorNumber == ROOT_PROCESSOR) printf("Checkpoint on CYCLE flag.\n");

#ifdef USE_HDF5_GROUPS

#ifdef RAD_HYDRO
      Group_WriteAllData(MetaData.HistoryDumpName, MetaData.CycleNumber, &TopGrid, MetaData, Exterior, ImplicitSolver);
#else
      Group_WriteAllData(MetaData.HistoryDumpName, MetaData.CycleNumber, &TopGrid, MetaData, Exterior);
#endif

#else

#ifdef RAD_HYDRO
      WriteAllData(MetaData.HistoryDumpName, MetaData.CycleNumber, &TopGrid, MetaData, Exterior, ImplicitSolver);
#else
      WriteAllData(MetaData.HistoryDumpName, MetaData.CycleNumber, &TopGrid, MetaData, Exterior);
#endif

#endif

    }
#endif
    
    // Write out the level timers

#ifdef TRACE_LEVELS
      fprintf(levelPtr, "%8"ISYM" %9"ISYM" %9"ISYM" %9"ISYM" %9"ISYM" %9"ISYM" %9"ISYM" %9"ISYM" %9"ISYM" %9"ISYM" %9"ISYM"\n", 
        MetaData.CycleNumber,
        level_count[0], level_count[1], level_count[2], level_count[3], level_count[4],
        level_count[5], level_count[6], level_count[7], level_count[8], level_count[9]);

      fprintf(levelPtr, "%8"ISYM" %9.3e %9.3e %9.3e %9.3e %9.3e %9.3e %9.3e %9.3e %9.3e %9.3e\n",
        MetaData.CycleNumber,
        level_timer[0], level_timer[1], level_timer[2], level_timer[3], level_timer[4],
        level_timer[5], level_timer[6], level_timer[7], level_timer[8], level_timer[9]);

      fprintf(levelPtr, "%8"ISYM" %9.3e %9.3e %9.3e %9.3e %9.3e %9.3e %9.3e %9.3e %9.3e %9.3e\n",
        MetaData.CycleNumber,
        rebuild_timer[0], rebuild_timer[1], rebuild_timer[2], rebuild_timer[3], rebuild_timer[4],
        rebuild_timer[5], rebuild_timer[6], rebuild_timer[7], rebuild_timer[8], rebuild_timer[9]);

      fflush(levelPtr);
#endif

#ifdef USE_MPI 
    MPI_Barrier(MPI_COMM_WORLD);
    tlev1 = MPI_Wtime();
#endif
 
    /* Rebuild the grids from level 0. */

#ifdef USE_MPI
    treb0 = MPI_Wtime();
#endif

#ifdef MEM_TRACE
    MemInUse = mused();
    fprintf(memtracePtr, "Pre loop rebuild %8"ISYM"  %16"ISYM" \n", MetaData.CycleNumber, MemInUse);
    fflush(memtracePtr);
#endif
 
    if (ProblemType != 25) {
      // CountGrids(MetaData, LevelArray);
      if (RebuildHierarchy(&MetaData, LevelArray, 0) == FAIL) {
	fprintf(stderr, "Error in RebuildHierarchy.\n");
	return FAIL;
      }
    }

#ifdef MEM_TRACE
    MemInUse = mused();
    fprintf(memtracePtr, "Post loop rebuild %8"ISYM"  %16"ISYM" \n", MetaData.CycleNumber, MemInUse);
    fflush(memtracePtr);
#endif

#ifdef USE_MPI
    treb1 = MPI_Wtime();
#endif
 
    /* Add time and check stopping criteria (steps #21 & #22)
       (note the topgrid is also keeping its own time but this statement will
       keep the two in synch). */
 
    MetaData.Time += dt;
    MetaData.CycleNumber++;
	
    if (MetaData.Time >= MetaData.StopTime) {
      if (MyProcessorNumber == ROOT_PROCESSOR)
	printf("Stopping on top grid time limit.\n");
      Stop = TRUE;
    }
    if (MetaData.CycleNumber >= MetaData.StopCycle) {
      if (MyProcessorNumber == ROOT_PROCESSOR)
	printf("Stopping on top grid cycle limit.\n");
      Stop = TRUE;
    }
 
    /* Check for time-actions. */
 
    if (CheckForTimeAction(LevelArray, MetaData) == FAIL) {
      fprintf(stderr, "Error in CheckForTimeActions.\n");
      return FAIL;
    }
 
    /* Check for output. */

#ifdef RAD_HYDRO 
    if (CheckForOutput(&TopGrid, MetaData, Exterior, ImplicitSolver, WroteData) == FAIL) {
      fprintf(stderr, "Error in CheckForOutput.\n");
      return FAIL;
    }
#else
    if (CheckForOutput(&TopGrid, MetaData, Exterior, WroteData) == FAIL) {
      fprintf(stderr, "Error in CheckForOutput.\n");
      return FAIL;
    }
#endif
 
    /* Try to cut down on memory fragmentation. */
 
#ifdef REDUCE_FRAGMENTATION
 
    if (WroteData && !Stop)
      if (ReduceFragmentation(TopGrid, MetaData, Exterior, LevelArray) == FAIL) {
	fprintf(stderr, "Error in ReduceFragmentation.\n");
	return FAIL;
      }
 
#endif /* REDUCE_FRAGMENTATION */

#ifdef USE_JBPERF
    jbPerf.stop("EL");
    if (((jb_iter+1) % JB_ITER_PER_SEGMENT)==0) jbPerf.end("EL");
#endif

#ifdef MEM_TRACE
    MemInUse = mused();
    fprintf(memtracePtr, "Bot %8"ISYM"  %16"ISYM" \n", MetaData.CycleNumber, MemInUse);
    fflush(memtracePtr);
#endif

  for ( i = 0; i < MAX_NUMBER_OF_TASKS; i++ ) {
    TaskMemory[i] = -1;
  }

#ifdef MEM_TRACE

  MPI_Datatype DataTypeInt = (sizeof(Eint64) == 4) ? MPI_INT : MPI_LONG_LONG_INT;
  MPI_Arg ThisTask;
  MPI_Arg TaskCount;
  MPI_Arg Count = 1;
  MPI_Arg stat;

  stat = MPI_Comm_size(MPI_COMM_WORLD, &TaskCount);
  stat = MPI_Comm_rank(MPI_COMM_WORLD, &ThisTask);
  stat = MPI_Allgather(&MemInUse, Count, DataTypeInt, TaskMemory, Count, DataTypeInt, MPI_COMM_WORLD);

/*
  if (ThisTask == 0 ) {
    for ( i = 0; i < TaskCount; i++) {
      fprintf(stderr, "TaskMemory : Task %"ISYM"  Memory %"ISYM"\n", i, TaskMemory[i]);
    }
  }
*/

#endif /* MEM_TRACE */

#ifdef USE_MPI
    tloop1 = MPI_Wtime();
#endif

    FILE *evlog;

    if (MyProcessorNumber == ROOT_PROCESSOR) {
      FLOAT z;
      z = 0.0;
      if (ComovingCoordinates) {
        FLOAT evt_a, evt_dadt;
        CosmologyComputeExpansionFactor(MetaData.Time, &evt_a, &evt_dadt);
        z = (1 + InitialRedshift)/evt_a - 1;
      }

      evlog = fopen("Evtime", "a");
      fprintf(evlog, "%8"ISYM"  %13.6e  %13.6e  %13.6e  %13.6e  %13.6e\n", MetaData.CycleNumber, MetaData.Time, z, tlev1-tlev0, treb1-treb0, tloop1-tloop0);
      fclose(evlog);
    }

#ifdef MEM_TRACE_END
    if (WroteData) {
      if (MemInUse > MemoryLimit) {
        if (MyProcessorNumber == ROOT_PROCESSOR)
          printf("Stopping due to memory limit.\n");
        Stop = TRUE;
      }
    }
#endif

    HaltFlag = CheckForHalt();
    if (HaltFlag) {
      if (MyProcessorNumber == ROOT_PROCESSOR)
        printf("Stopping on HALT flag.\n");
      Stop = TRUE;
    }
 
  } // ===== end of main loop ====
 
#ifdef USE_JBPERF
  if (((jb_iter+1) % JB_ITER_PER_SEGMENT)!=0) jbPerf.end("EL");
  jbPerf.attribute ("timestep",0, JB_NULL);
#endif

#ifdef USE_MPI
  MetaData.CPUTime = MPI_Wtime() - MetaData.CPUTime;
#endif
 
  /* Done, so report on current time, etc. */
 
  if (MyProcessorNumber == ROOT_PROCESSOR) {
    printf("Time     = %9"FSYM"   CycleNumber = %6"ISYM"    Wallclock   = %9"FSYM"\n",
	   MetaData.Time, MetaData.CycleNumber, MetaData.CPUTime);
    printf("StopTime = %9"FSYM"   StopCycle   = %6"ISYM"\n",
	   MetaData.StopTime, MetaData.StopCycle);
  }
 
  /* If we are running problem 23, TestGravity, then check the results. */
 
  if (ProblemType == 23)
    TestGravityCheckResults(LevelArray);
  if (ProblemType == 25 && NumberOfProcessors == 0)
    TestGravitySphereCheckResults(LevelArray);
 
  /* if we are doing data dumps, then do one last one */
 
  if (MetaData.dtDataDump != 0.0 || MetaData.CycleSkipDataDump != 0 || HaltFlag != 0) {

#ifdef USE_HDF5_GROUPS

#ifdef RAD_HYDRO
      Group_WriteAllData(MetaData.DataDumpName, MetaData.DataDumpNumber, &TopGrid, MetaData, Exterior, ImplicitSolver, -666);
#else
      Group_WriteAllData(MetaData.DataDumpName, MetaData.DataDumpNumber, &TopGrid, MetaData, Exterior, -666);
#endif 
      
#else

#ifdef RAD_HYDRO
      WriteAllData(MetaData.DataDumpName, MetaData.DataDumpNumber, &TopGrid, MetaData, Exterior, ImplicitSolver, -666);
#else
      WriteAllData(MetaData.DataDumpName, MetaData.DataDumpNumber, &TopGrid, MetaData, Exterior, -666);
#endif

#endif

  }

  if (NumberOfProcessors > 1)
    printf("Communication: processor %"ISYM" CommunicationTime = %"FSYM"\n",
	   MyProcessorNumber, CommunicationTime);
 
  /* done */

#ifdef USE_MPI
  texit = MPI_Wtime();
#endif
 
  return SUCCESS;
}
