/***********************************************************************
/
/  COMMUNICATION ROUTINE: LOAD BALANCE GRIDS
/
/  written by: Greg Bryan
/  date:       December, 1997
/  modified1:
/
/  PURPOSE:
/
************************************************************************/
 
#include <stdio.h>
#include <string.h>
#ifdef USE_MPI
#include "mpi.h"
#endif
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "TopGridData.h"
#include "Hierarchy.h"
#include "LevelHierarchy.h"
 
// Function prototypes
 
void WriteListOfFloats(FILE *fptr, int N, float floats[]);
void fpcol(float *x, int n, int m, FILE *fptr);
 
#define LOAD_BALANCE_RATIO 1.05
 
int CommunicationLoadBalanceGrids(HierarchyEntry *GridHierarchyPointer[],
				  int NumberOfGrids)
{
 
  if (NumberOfProcessors == 1 || NumberOfGrids <= 1)
    return SUCCESS;
 
  /* Initialize */
 
  int i, GridMemory, NumberOfCells, CellsTotal, Particles;
  float AxialRatio, GridVolume;
  float *ComputeTime = new float[NumberOfGrids];
  float *ProcessorComputeTime = new float[NumberOfProcessors];
 
#ifdef MPI_INSTRUMENTATION
  starttime = MPI_Wtime();
  moving_count ++;
  out_count ++;
#endif
 
  for (i = 0; i < NumberOfProcessors; i++)
    ProcessorComputeTime[i] = 0;
 
  /* Compute work for each grid. */
 
  for (i = 0; i < NumberOfGrids; i++) {
    GridHierarchyPointer[i]->GridData->CollectGridInformation
      (GridMemory, GridVolume, NumberOfCells, AxialRatio, CellsTotal, Particles);
    ComputeTime[i] = GridMemory; // roughly speaking
    ProcessorComputeTime[GridHierarchyPointer[i]->GridData->ReturnProcessorNumber()] += GridMemory;
  }
 
  /* Transfer grids from heavily-loaded processors. */
 
  int Done = FALSE, MinProc = 0, MaxProc = 0;
  while (!Done) {
 
    /* Find min and max */
 
    float MaxVal = 0, MinVal = huge_number;

    MaxProc = -1;
    MinProc = -1;

    //dcc 09/22/05 updated this loop to avoid huge_number being too small.

    for (i = 0; i < NumberOfProcessors; i++) {
      if (ProcessorComputeTime[i] > MaxVal) {
	MaxVal = ProcessorComputeTime[i];
	MaxProc = i;
      }
    }
    for (i = 0; i < NumberOfProcessors; i++) {
      if (ProcessorComputeTime[i] < MinVal) {
	MinVal = ProcessorComputeTime[i];
	MinProc = i;
      }
    }

    if(MinProc == -1 || MaxProc == -1 )
      fprintf(stderr, "TERRIBLE ERROR: CommunicationLoadBalance unable to find processors.\n");


    /* Transfer a grid if the ratio is large enough. */
 
    if (MaxVal > LOAD_BALANCE_RATIO*MinVal) {
 
      /* Find a grid to transfer. */
 
      for (i = 0; i < NumberOfGrids; i++) {
	int proc = GridHierarchyPointer[i]->GridData->ReturnProcessorNumber();
	if (proc == MaxProc && ComputeTime[i] < 0.5*(MaxVal-MinVal)) {
 
	  /* Transfer. */
 
//	  printf("%"ISYM": moving grid %"ISYM" from %"ISYM" -> %"ISYM"\n", MyProcessorNumber, i, proc, MinProc);
 
          /* Attach ForcingFields before transfer, if necessary; then detach */
 
          if (RandomForcing)  //AK
            GridHierarchyPointer[i]->GridData->AppendForcingToBaryonFields();
          GridHierarchyPointer[i]->GridData->CommunicationMoveGrid(MinProc);
          if (RandomForcing)  //AK
            GridHierarchyPointer[i]->GridData->RemoveForcingFromBaryonFields();
 
//	  printf("%"ISYM": done moving grid %"ISYM"\n", MyProcessorNumber, i);
 
	  /* Update processor compute times. */
 
	  ProcessorComputeTime[MaxProc] -= ComputeTime[i];
	  ProcessorComputeTime[MinProc] += ComputeTime[i];
 
	  break;
	}
      }
 
      /* If we didn't find an appropriate transfer then quit. */
 
      if (i == NumberOfGrids) {
	Done = TRUE;
#ifdef MPI_INSTRUMENTATION
	if (MinVal == 0)
	  timer[3] = MaxVal;
	else
	  timer[3] = MaxVal/MinVal;
#endif /* MPI_INSTRUMENTATION */
      }
    }
    else {
      Done = TRUE;
#ifdef MPI_INSTRUMENTATION
      counter[3]++;
      if (MinVal == 0)
	timer[3] = MaxVal;
      else
	timer[3] = MaxVal/MinVal;
#endif /* MPI_INSTRUMENTATION */
    }
  }
 
#ifdef MPI_INSTRUMENTATION
  moving_pct += float(out_count)/NumberOfGrids;
#endif /* MPI_INSTRUMENTATION */
 
  if (MyProcessorNumber == ROOT_PROCESSOR) {
    printf("LoadBalance (grids=%"ISYM"): \n", NumberOfGrids);
    float norm = ProcessorComputeTime[0];
    for (i = 1; i < NumberOfProcessors; i++)
      norm = max(norm, ProcessorComputeTime[i]);
    for (i = 0; i < NumberOfProcessors; i++)
      ProcessorComputeTime[i] /= max(norm, 1.0e-10);
    // WriteListOfFloats(stdout, NumberOfProcessors, ProcessorComputeTime);
    fpcol(ProcessorComputeTime, NumberOfProcessors, 16, stdout);
  }
 
  delete [] ComputeTime;
  delete [] ProcessorComputeTime;
 
#ifdef MPI_INSTRUMENTATION
  endtime = MPI_Wtime();
  timer[2] += endtime - starttime;
  counter[2] ++;
#endif /* MPI_INSTRUMENTATION */
 
  return SUCCESS;
}
