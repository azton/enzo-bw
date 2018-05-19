/***********************************************************************
/
/  AMR MAIN CODE
/
/  written by: Greg Bryan
/  date:       November, 1994
/  modified:   Robert Harkness
/  date:       August 12th 2006
/              April 27th 2009
/
/  PURPOSE:
/    This is main() for the amr code.  It interprets the arguments and
/    then calls the appropriate routines depending on whether we are
/    doing a new run, a restart, an extraction, or a projection.
/
************************************************************************/

#ifdef RAD_HYDRO
#include "ImplicitProblemABC_preincludes.h"
#endif
 
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
 
#ifdef USE_MPI
#include "mpi.h"
#endif /* USE_MPI */
 
#include "version.def"
#include "performance.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#define DEFINE_STORAGE
#include "global_data.h"
#include "units.h"
#include "flowdefs.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "Hierarchy.h"
#include "LevelHierarchy.h"
#include "TopGridData.h"
#include "CosmologyParameters.h"
#include "StarParticleData.h"
#ifdef RAD_HYDRO
#include "ImplicitProblemABC.h"
#include "gFLDProblem.h"
#include "gFLDSplit.h"
#include "DualFLD.h"
#include "FSProb.h"
#include "NullProblem.h"
#endif
#undef DEFINE_STORAGE
 
// Function prototypes
 
int InitializeNew(  char *filename, HierarchyEntry &TopGrid, TopGridData &tgd,
		    ExternalBoundary &Exterior, float *Initialdt);

int InitializeLocal(int restart, HierarchyEntry &TopGrid, 
		    TopGridData &MetaData
#ifdef RAD_HYDRO
		    , ImplicitProblemABC *ImplicitSolver
#endif
		    );


int ReadAllData(char *filename, HierarchyEntry *TopGrid, TopGridData &tgd,
		    ExternalBoundary *Exterior);
int Group_ReadAllData(char *filename, HierarchyEntry *TopGrid, TopGridData &tgd,
		    ExternalBoundary *Exterior);

int EvolveHierarchy(HierarchyEntry &TopGrid, TopGridData &tgd,
		    ExternalBoundary *Exterior, 
#ifdef RAD_HYDRO
		    ImplicitProblemABC *ImplicitSolver,
#endif
		    LevelHierarchyEntry *Array[], float Initialdt);

void ExtractSection(HierarchyEntry &TopGrid, TopGridData &tgd,
		    LevelHierarchyEntry *Array[], ExternalBoundary *Exterior,
		    int ExtractStart[], int ExtractEnd[],
		    FLOAT ProjectStartCoordinates[],
		    FLOAT ProjectEndCoordinates[], int ExtractLevel);
int OutputLevelInformation(FILE *fptr, TopGridData &tgd,
			   LevelHierarchyEntry *Array[]);
int ProjectToPlane(TopGridData &MetaData, LevelHierarchyEntry *LevelArray[],
		   int ProjectStart[], int ProjectEnd[],
		   FLOAT ProjectStartCoordinates[],
		   FLOAT ProjectEndCoordinates[], int ProjectLevel,
		   int ProjectionDimension, char *ProjectionFileName,
		   int ProjectionSmooth, ExternalBoundary *Exterior);
int OutputAsParticleData(TopGridData &MetaData,
			 LevelHierarchyEntry *LevelArray[],
			 int RegionStart[], int RegionEnd[],
			 FLOAT RegionStartCoordinates[],
			 FLOAT RegionEndCoordinates[], int RegionLevel,
			 char *OutputFileName);
int InterpretCommandLine(int argc, char *argv[], char *myname,
			 int &restart, int &debug, int &extract,
			 int &InformationOutput,
			 int &OutputAsParticleDataFlag,
			 int &project, int &ProjectionDimension,
			 int &ProjectionSmooth,
			 char *ParameterFile[],
			 int RegionStart[], int RegionEnd[],
			 FLOAT RegionStartCoordinates[],
			 FLOAT RegionEndCoordinates[],
			 int &Level, int MyProcessorNumber);
void AddLevel(LevelHierarchyEntry *Array[], HierarchyEntry *Grid, int level);
int SetDefaultGlobalValues(TopGridData &MetaData);
int CommunicationInitialize(Eint32 *argc, char **argv[]);
int CommunicationFinalize();
void CommunicationAbort(int);
int CommunicationPartitionGrid(HierarchyEntry *Grid);
int ENZO_OptionsinEffect(void);

#ifdef TASKMAP
int GetNodeFreeMemory(void);
#endif

void my_exit(int status);

#ifdef RAD_HYDRO
int DetermineParallelism(HierarchyEntry *TopGrid, TopGridData &MetaData);
#endif
 
#ifdef MEM_TRACE
Eint64 mused(void);
#endif 


 
//  ENZO Main Program

 
Eint32 main(Eint32 argc, char *argv[])
{

#ifdef MEM_TRACE
    Eint64 MemInUse;
#endif

  int i;

  // Initialize Communications

  CommunicationInitialize(&argc, &argv);

  int int_argc;
  int_argc = argc;
 
  if (MyProcessorNumber == ROOT_PROCESSOR) {
    printf("=========================\n");
    printf("ENZO CVS VERSION %d.%d.%d\n",	   
	   ENZO_VERSION_MAJOR,
	   ENZO_VERSION_MINOR,
	   ENZO_VERSION_REVISION);
    printf("=========================\n");
    fflush(stdout);
  }


  // Performance Monitoring

#ifdef USE_MPI
  double t_init0, t_init1;

  t_init0 = MPI_Wtime();
#endif /* USE_MPI */

  // Task Mapping

  for (i = 0; i < MAX_NUMBER_OF_TASKS; i++ ) {
    TaskMemory[i] = -1;
  }

#ifdef TASKMAP
  GetNodeFreeMemory();
#endif


  // Begin 

#ifdef MPI_INSTRUMENTATION
  Start_Wall_Time = MPI_Wtime();
#endif

#ifdef OOC_BOUNDARY
  ExternalBoundaryIO = TRUE;
  ExternalBoundaryTypeIO = FALSE;
  ExternalBoundaryValueIO = FALSE;
#else
  ExternalBoundaryIO = FALSE;
  ExternalBoundaryTypeIO = FALSE;
  ExternalBoundaryValueIO = FALSE;
#endif

  ENZO_OptionsinEffect();


#ifdef ISOLATED_GRAVITY
#ifdef UNIGRID_TRANSPOSE
    /* it turns out that Robert Harkness' unigrid transpose stuff is incompatible with the top
       grid isolated gravity boundary conditions.  I'm not 100 percent sure why this is - in the 
       meantime, just double-check to make sure that if one tries to use the isolated boundary
       conditions when the unigrid transpose stuff is on, the code crashes loudly.
       -- BWO, 26 June 2008 */
      if (MyProcessorNumber == ROOT_PROCESSOR){
	fprintf(stderr, "\n\n");
	fprintf(stderr, "  ************************************************************************\n");
	fprintf(stderr, "  ****  D'oh!  At present, you cannot use isolated top grid boundary  ****\n");
	fprintf(stderr, "  ****  conditions with the top grid unigrid bookkeeping scheme.      ****\n");
	fprintf(stderr, "  ****  Consult Brian O'Shea for the details of this wackiness,       ****\n");
	fprintf(stderr, "  ****  and in the meantime recompile enzo with                       ****\n");
	fprintf(stderr, "  ****  'gmake unigrid-transpose-no' instead of                       ****\n");
	fprintf(stderr, "  ****  'gmake unigrid-transpose-yes'.                                ****\n");
	fprintf(stderr, "  ************************************************************************\n");      
	fprintf(stderr, "\n\n");
      }
      my_exit(EXIT_FAILURE);
#endif /* UNIGRID_TRANSPOSE */
#endif /* ISOLATED_GRAVITY */



  // Main declarations
 
  TopGridData MetaData;
  HierarchyEntry TopGrid;
  ExternalBoundary Exterior;
  LevelHierarchyEntry *LevelArray[MAX_DEPTH_OF_HIERARCHY];

 
  // Initialize
 
  int restart                  = FALSE,
      OutputAsParticleDataFlag = FALSE,
      InformationOutput        = FALSE,
      project                  = FALSE,
      ProjectionDimension      = INT_UNDEFINED,
      ProjectionSmooth         = FALSE;
  debug                        = FALSE;
  extract                      = FALSE;
  flow_trace_on                = FALSE;
  char *ParameterFile          = NULL;
  char *myname                 = argv[0];

  int RegionStart[MAX_DIMENSION],
      RegionEnd[MAX_DIMENSION],
      RegionLevel;
  FLOAT RegionStartCoordinates[MAX_DIMENSION],
        RegionEndCoordinates[MAX_DIMENSION];

  float Initialdt              = 0;
 
#ifdef FLOW_TRACE
  char pid[MAX_TASK_TAG_SIZE];
  char flow_trace_log[MAX_NAME_LENGTH];
  sprintf(pid, "%"TASK_TAG_FORMAT""ISYM, MyProcessorNumber);
  strcpy(flow_trace_log, "FlowTrace.");
  strcat(flow_trace_log, pid);
  flow_trace_fptr = fopen( flow_trace_log, "w" );
  flow_trace_level = 0;
  flow_trace_on = TRUE;
#endif
 
  if (flow_trace_on) flow_trace1("ENZO");
 
#ifdef MPI_INSTRUMENTATION
  char perfname[MAX_NAME_LENGTH];
 
  flagging_count = 0;
  moving_count = 0;
  in_count = 0;
  out_count = 0;
  flagging_pct = 0.0;
  moving_pct = 0.0;
 
  GlobalCommunication = 0.0;
  RecvComm = 0.0;
  WaitComm = 0.0;
 
  for (i=0; i < MAX_COUNTERS; i++) {
    timer[i] = 0.0;
    counter[i] = 0;
  }
 
  sprintf(perfname, "perfdata_%"ISYM".%"ISYM, NumberOfProcessors, MyProcessorNumber);
  perfname[strlen(perfname)] = '\0';

  filePtr = fopen(perfname, "w");
 
  sprintf(tracename, "trace_%"ISYM".%"ISYM, NumberOfProcessors, MyProcessorNumber);
  tracename[strlen(tracename)] = '\0';

  sprintf(level_timer_name, "levels_%"ISYM".%"ISYM, NumberOfProcessors, MyProcessorNumber);
  level_timer_name[strlen(level_timer_name)] = '\0';

#ifdef TRACE_LEVELS
  levelPtr = fopen(level_timer_name, "w");
#endif

#ifdef MPI_TRACE
  tracePtr = fopen(tracename, "w");
  traceMPI = TRUE;
#else
  traceMPI = FALSE;
#endif /* MPI_TRACE */

#endif /* MPI_INSTRUMENTATION */

#ifdef MEM_TRACE
  sprintf(memtracename, "mem_%"ISYM".%"ISYM, NumberOfProcessors, MyProcessorNumber);
  memtracename[strlen(memtracename)] = '\0';
  memtracePtr = fopen(memtracename, "a");
  traceMEM = TRUE;
#else
  traceMEM = FALSE;
#endif /* MEM_TRACE */




  // START
 
  for (int level = 0; level < MAX_DEPTH_OF_HIERARCHY; level++) {
    LevelArray[level] = NULL;
  }
 

  // Interpret command-line arguments
 
  if (InterpretCommandLine(int_argc, argv, myname, restart, debug, extract,
			   InformationOutput, OutputAsParticleDataFlag,
			   project, ProjectionDimension, ProjectionSmooth,
			   &ParameterFile,
			   RegionStart, RegionEnd,
			   RegionStartCoordinates, RegionEndCoordinates,
			   RegionLevel, MyProcessorNumber) == FAIL) {
    if(int_argc==1){
      my_exit(EXIT_SUCCESS);
    } else {
      my_exit(EXIT_FAILURE);
    }
  }

 
  // If we need to read the parameter file as a restart file, do it now
 
  if (restart || OutputAsParticleDataFlag || extract || InformationOutput || project) {
 
    SetDefaultGlobalValues(MetaData);
 
    if (debug) printf("Reading parameter file %s\n", ParameterFile);

#ifdef USE_HDF5_GROUPS
    if (debug) fprintf(stderr, "Input with Group_ReadAllData\n");
    if (Group_ReadAllData(ParameterFile, &TopGrid, MetaData, &Exterior) == FAIL) {
      if (MyProcessorNumber == ROOT_PROCESSOR)
	fprintf(stderr, "Error in ParameterFile %s.\n", ParameterFile);
      my_exit(EXIT_FAILURE);
    }
#else 
    if (debug) fprintf(stderr, "Input with ReadAllData\n");
    if (ReadAllData(ParameterFile, &TopGrid, MetaData, &Exterior) == FAIL) {
      if (MyProcessorNumber == ROOT_PROCESSOR)
	fprintf(stderr, "Error in ParameterFile %s.\n", ParameterFile);
      my_exit(EXIT_FAILURE);
    }
#endif
 
    if (!ParallelRootGridIO && restart && TopGrid.NextGridThisLevel == NULL) {
      CommunicationPartitionGrid(&TopGrid);  // partition top grid if necessary
    }
 
    if (MyProcessorNumber == ROOT_PROCESSOR) {
      fprintf(stderr, "Successfully read ParameterFile %s.\n", ParameterFile);
    }
 
    AddLevel(LevelArray, &TopGrid, 0);    // recursively add levels
 
  }


  // Information dump
 
  if (InformationOutput) {
    OutputLevelInformation(stdout, MetaData, LevelArray);
    my_exit(EXIT_SUCCESS);
  }
 
 
  // Extract a grid subsection (doesn't return)
 
  if (extract) {
    ExtractSection(TopGrid, MetaData, LevelArray, &Exterior,
		   RegionStart, RegionEnd,
		   RegionStartCoordinates, RegionEndCoordinates, RegionLevel);
  }
 
 
  // Output fields as particle data
 
  if (OutputAsParticleDataFlag) {
    if (OutputAsParticleData(MetaData, LevelArray, RegionStart, RegionEnd,
			     RegionStartCoordinates, RegionEndCoordinates,
			     RegionLevel, "amr.particles") == FAIL)
      my_exit(EXIT_FAILURE);
    else
      my_exit(EXIT_SUCCESS);
  } 

 
  // Project 3D field to 2D plane
 
  if (project) {
    if (ProjectToPlane(MetaData, LevelArray, RegionStart, RegionEnd,
		     RegionStartCoordinates, RegionEndCoordinates,
		     RegionLevel, ProjectionDimension, "amr.project",
		     ProjectionSmooth, &Exterior) == FAIL)
      my_exit(EXIT_FAILURE);
    else
      my_exit(EXIT_SUCCESS);
  } 

 
  // Normal start: Open and read parameter file

#ifdef MEM_TRACE
    MemInUse = mused();
    fprintf(memtracePtr, "Call initialize %16"ISYM" \n", MemInUse);
#endif

 
  if (!restart) {

    if (InitializeNew(ParameterFile, TopGrid, MetaData, Exterior, &Initialdt) == FAIL) {
      if (MyProcessorNumber == ROOT_PROCESSOR)
	fprintf(stderr, "Error in Parameter File %s.\n", ParameterFile);
      my_exit(EXIT_FAILURE);
    }
    else {
      if (MyProcessorNumber == ROOT_PROCESSOR)
	fprintf(stderr, "Successfully read in parameter file %s.\n", ParameterFile);
      AddLevel(LevelArray, &TopGrid, 0);    // recursively add levels
    }
 
#ifdef USE_MPI
    MPI_Barrier(MPI_COMM_WORLD);
    t_init1 = MPI_Wtime();
    fprintf(stderr, "INITIALIZATION TIME = %16.8e\n", (t_init1-t_init0));
    MPI_Barrier(MPI_COMM_WORLD);
#endif /* USE_MPI */

  }

#ifdef RAD_HYDRO
  // Determine Top-Grid Parallelism information (store in grids)
  if (DetermineParallelism(&TopGrid, MetaData) == FAIL) {
    fprintf(stderr,"Error in DetermineParallelism.\n");
    my_exit(EXIT_FAILURE);
  }
#endif

  // if using an implicit RT solver, declare the appropriate object here
#ifdef RAD_HYDRO
  ImplicitProblemABC *ImplicitSolver;
  if (ImplicitProblem == 1) {
    ImplicitSolver = new gFLDProblem; }
  else if (ImplicitProblem == 2) {
    ImplicitSolver = new FSProb; }
  else if (ImplicitProblem == 3) {
    ImplicitSolver = new gFLDSplit; }
  else if (ImplicitProblem == 4) {
    ImplicitSolver = new DualFLD; }
  else {
    ImplicitSolver = new NullProblem; }
#endif

  // Perform local initialization (even for restart, may just return, depends on problem)

  if (InitializeLocal(restart, TopGrid, MetaData
#ifdef RAD_HYDRO
		     , ImplicitSolver
#endif
		     ) == FAIL) {
    if (MyProcessorNumber == ROOT_PROCESSOR)
      fprintf(stderr, "Error in Local Initialization.\n");
    my_exit(EXIT_FAILURE); 
  }


#ifdef MEM_TRACE
  MemInUse = mused();
  fprintf(memtracePtr, "Call evolve hierarchy %8"ISYM"  %16"ISYM" \n", MetaData.CycleNumber, MemInUse);
#endif



 
  // Call the main evolution routine
 
  if (EvolveHierarchy(TopGrid, MetaData, &Exterior,
#ifdef RAD_HYDRO
		      ImplicitSolver,
#endif
		      LevelArray, Initialdt) == FAIL) {
    if (MyProcessorNumber == ROOT_PROCESSOR) {
      fprintf(stderr, "Error in EvolveHierarchy.\n");
    }
    my_exit(EXIT_FAILURE);
  }
  else
  {
    if (MyProcessorNumber == ROOT_PROCESSOR) {
      fprintf(stderr, "Successful run, exiting.\n");
    }
  }

 
  if (flow_trace_on) flow_trace2("ENZO");




  // Terminate
  
#ifdef MPI_INSTRUMENTATION

  End_Wall_Time = MPI_Wtime();
  WallTime = End_Wall_Time - Start_Wall_Time;
 
  fprintf(filePtr, "Elapsed wall time:                   %12.6e\n", WallTime);
  fprintf(filePtr, "Communication time:                  %12.6e\n", CommunicationTime);
  fprintf(filePtr, "Global communication time:           %12.6e\n", GlobalCommunication);
  fprintf(filePtr, "Receive communication time:          %12.6e\n", RecvComm);
  fprintf(filePtr, "Waiting communication time:          %12.6e\n", WaitComm);
 
  fprintf(filePtr, "\n\n");
  fprintf(filePtr, "Transferring region       (%8"ISYM" times) %12.6e\n", counter[5],  timer[5]);
  fprintf(filePtr, "Sending particles         (%8"ISYM" times) %12.6e\n", counter[7],  timer[7]);
  fprintf(filePtr, "Transferring particles    (%8"ISYM" times) %12.6e\n", counter[9],  timer[9]);
  fprintf(filePtr, "Transferring Fluxes       (%8"ISYM" times) %12.6e\n", counter[12], timer[12]);
  fprintf(filePtr, "ShareGrids                (%8"ISYM" times) %12.6e\n", counter[13], timer[13]);
  fprintf(filePtr, "Transpose                 (%8"ISYM" times) %12.6e\n", counter[14], timer[14]);
  fprintf(filePtr, "BroadcastValue            (%8"ISYM" times) %12.6e\n", counter[15], timer[15]);
  fprintf(filePtr, "MinValue                  (%8"ISYM" times) %12.6e\n", counter[16], timer[16]);
  fprintf(filePtr, "UpdateStarParticleCount   (%8"ISYM" times) %12.6e\n", counter[11], timer[11]);
 
  fprintf(filePtr, "\n\n");
  fprintf(filePtr, "RebuildHierarchy          (%8"ISYM" times) %12.6e\n", counter[1],  timer[1]);
  fprintf(filePtr, "RebuildHierarchy interval (%8"ISYM" times) %12.6e\n", counter[0],  timer[0]);
  fprintf(filePtr, "Load balancing            (%8"ISYM" times) %12.6e\n", counter[2],  timer[2]);
  fprintf(filePtr, "Region transfer size      (%8"ISYM" times) %12.6e\n", counter[5],  timer[6]);
  fprintf(filePtr, "Particles sent            (%8"ISYM" times) %12.6e\n", counter[7],  timer[8]);
  fprintf(filePtr, "Particle transfer size    (%8"ISYM" times) %12.6e\n", counter[9],  timer[10]);
 
  fprintf(filePtr, "\n\n");
  fprintf(filePtr, "Number of load balancing calls %"ISYM"/%"ISYM" (LOAD_BALANCE_RATIO=%"FSYM")\n",counter[3], counter[2], timer[3]);
  fprintf(filePtr, "Number of flagging cells  (%8"ISYM" times) %12.6e\n", counter[4], timer[4]);
 
  fprintf(filePtr, "\n\n");

  if ( flagging_count != 0 )
    fprintf(filePtr, "Average percentage of flagging cells %12.6e(= %12.6e/%"ISYM")\n",
            flagging_pct/flagging_count, flagging_pct, flagging_count);
  else
    fprintf(filePtr, "Average percentage of flagging cells 0\n");
 
  if ( moving_count != 0 )
    fprintf(filePtr, "Average percentage of moving cells   %12.6e(= %12.6e/%"ISYM")\n",
            moving_pct/moving_count, moving_pct,moving_count);
  else
    fprintf(filePtr, "Average percentage of moving cells 0\n");
 
  fclose(filePtr);

#ifdef TRACE_LEVELS
  fclose(levelPtr);
#endif
 
#ifdef MPI_TRACE
  fclose(tracePtr);
#endif /* MPI_TRACE */

#endif /* MPI_INSTRUMENTATION */

#ifdef MEM_TRACE
  fclose(memtracePtr);
#endif /* MEM_TRACE */
 
#ifdef FLOW_TRACE
    fclose(flow_trace_fptr);
#endif
 
  my_exit(EXIT_SUCCESS);
 
}
 
 
 
 
void my_exit(int status)
{
  // Exit gracefully if successful; abort on error

  if (status == EXIT_SUCCESS) {

    if (MyProcessorNumber==0) {
      fprintf (stdout,"%s:%d Exiting.\n", __FILE__,__LINE__);
    }

    CommunicationFinalize();

    exit(status);

  } else if (status == EXIT_FAILURE) {

    fprintf (stderr,"%s:%d %"ISYM" ABORT ON EXIT_FAILURE!\n",
	     __FILE__,__LINE__,MyProcessorNumber);

    CommunicationAbort(status);

  } else {

    fprintf (stderr,"%s:%d %"ISYM" ABORT ON UNKNOWN EXIT VALUE %"ISYM"!\n",
	     __FILE__,__LINE__,MyProcessorNumber,status);

    CommunicationAbort(status);

  }
}
