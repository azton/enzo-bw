/***********************************************************************
/
/  COMMUNICATION ROUTINE: PARTITION GRID
/
/  written by: Greg Bryan
/  date:       December, 1997
/  modified1:  Robert Harkness
/  date:       March, 2006
/
/  PURPOSE:
/
************************************************************************/
 
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
 
#ifdef USE_MPI
#include "mpi.h"
#endif /* USE_MPI */
 
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
 
int CommunicationBroadcastValue(int *Value, int BroadcastProcessor);
int Enzo_Dims_create(int nnodes, int ndims, int *dims);
 
#define USE_OLD_CPU_DISTRIBUTION
 
 
int CommunicationPartitionGrid(HierarchyEntry *Grid)
{
 
  if (NumberOfProcessors == 1)
    return SUCCESS;
 
  // Declarations
 
  int Rank, dim, i, j, k, ijk, Dims[MAX_DIMENSION], Layout[] = {0,0,0};
  int TempDims[MAX_DIMENSION] /* , TempStart[MAX_DIMENSION] */ ;
  int *GridDims[MAX_DIMENSION], *StartIndex[MAX_DIMENSION];
  FLOAT Left[MAX_DIMENSION], Right[MAX_DIMENSION];
  FLOAT LeftEdge[MAX_DIMENSION], RightEdge[MAX_DIMENSION];
 
  //  printf("Enter CommunicationPartitionGrid\n");

#ifdef USE_MPI
  MPI_Barrier(MPI_COMM_WORLD);
#endif /* USE_MPI */
 
  /* Attach RandomForcingFields as BaryonFields (for the duration
     of partitioning only). */
 
  if (RandomForcing == 1 && ParallelRootGridIO != 1)
    Grid->GridData->AppendForcingToBaryonFields(); //AK
 
  /* Compute side length of root grid. */
 
  Grid->GridData->ReturnGridInfo(&Rank, Dims, Left, Right);
 
  for (dim = 0; dim < Rank; dim++)
    Dims[dim] -= 2*DEFAULT_GHOST_ZONES;
 
  float Edge = POW(float(Dims[0]*Dims[1]*Dims[2])/float(NumberOfProcessors),
		   1/float(Rank));
 
 
  /* If using MPI, use their routine to calculate layout. */
 
#ifdef USE_MPI
 
  int LayoutTemp[] = {0,0,0};

/*
  MPI_Arg Nnodes = NumberOfProcessors;
  MPI_Arg Ndims = Rank;
  MPI_Arg LayoutDims[] = {0, 0, 0};
 
  if (MPI_Dims_create(Nnodes, Ndims, LayoutDims) != MPI_SUCCESS) {
    fprintf(stderr, "Error in MPI_Dims_create.\n");
    return FAIL;
  }
*/

  int Nnodes = NumberOfProcessors;
  int Ndims = Rank;
  int LayoutDims[] = {0, 0, 0};

  if (Enzo_Dims_create(Nnodes, Ndims, LayoutDims) != SUCCESS) {
    fprintf(stderr, "Error in Enzo_Dims_create.\n");
    return FAIL;
  }

  for (dim = 0; dim < Rank; dim++)
    LayoutTemp[dim] = LayoutDims[dim];
 
  /* Swap layout because we want smallest value to be at Layout[0]. */
 
  for (dim = 0; dim < Rank; dim++)
    Layout[dim] = LayoutTemp[Rank-1-dim];
 
  /* Force some distributions if the default is brain-dead. */

/*
  if (Rank == 3 && NumberOfProcessors == 8)
    for (dim = 0; dim < Rank; dim++)
      Layout[dim] = 2;

  if (Rank == 3 && NumberOfProcessors == 64)
    for (dim = 0; dim < Rank; dim++)
      Layout[dim] = 4;

  if (Rank == 3 && NumberOfProcessors == 125)
    for (dim = 0; dim < Rank; dim++)
      Layout[dim] = 5;
 
  if (Rank == 3 && NumberOfProcessors == 216)
    for (dim = 0; dim < Rank; dim++)
      Layout[dim] = 6;
*/

  if (MyProcessorNumber == ROOT_PROCESSOR) {
    fprintf(stderr, "ENZO_layout %"ISYM" x %"ISYM" x %"ISYM"\n", Layout[0], Layout[1], Layout[2]);
  }

#endif /* USE_MPI */
 
 
  /* Generate arrays of grid dimensions and start positions. */
 
  int NumberOfNewGrids = 1;
 
  for (dim = 0; dim < MAX_DIMENSION; dim++) {
 
    /* Compute number of new grids along this dimension. */
 
    if (Layout[dim] == 0)
      Layout[dim] = max(nint((float)(Dims[dim])/Edge), 1);
 
    GridDims[dim] = new int[Layout[dim]];
    StartIndex[dim] = new int[Layout[dim]];
 
    /* Compute dims and start indexes of the grids along this dim. */
 
    float ExactDims = float(Dims[dim])/float(Layout[dim]);
    float ExactCount = 0;
    int DisplacementCount = 0;
 
    for (i = 0; i < Layout[dim]; i++) {
      ExactCount += ExactDims;
 
      /* Compute the integer number of cells along this dimension
	 (if dim == 0 then make sure it is even as well since the FFT
	 requires this). */
 
      if (dim == 0)
	GridDims[dim][i] = nint(ExactCount*0.5)*2 - DisplacementCount;
      else
	GridDims[dim][i] = nint(ExactCount) - DisplacementCount;
 
      StartIndex[dim][i] = DisplacementCount;
      DisplacementCount += GridDims[dim][i];
    }
 
    NumberOfNewGrids *= Layout[dim];
 
  }
 
  if (MyProcessorNumber == ROOT_PROCESSOR)
  {
    printf("PartitionGrid (on all processors): Layout = %"ISYM" %"ISYM" %"ISYM"\n",
      Layout[0], Layout[1], Layout[2]);
    printf("NumberOfNewGrids = %"ISYM"\n",NumberOfNewGrids);
 
    for (dim = 0; dim < MAX_DIMENSION; dim++)
    {
      printf("GridDims[%"ISYM"]: ",dim);
      for (i = 0; i < Layout[dim]; i++)
      {
        printf(" %"ISYM,GridDims[dim][i]);
      }
      printf("\n");
    }
    for (dim = 0; dim < MAX_DIMENSION; dim++)
    {
      printf("StartIndex[%"ISYM"]: ",dim);
      for (i = 0; i < Layout[dim]; i++)
      {
        printf(" %"ISYM,StartIndex[dim][i]);
      }
      printf("\n");
    }
  }
 
/*
  if ((ProblemType == 30) && (ParallelRootGridIO == 1) && (ParallelParticleIO == 1))
  {
    printf("Unigrid: %"ISYM"\n", Unigrid);
    printf("Set Unigrid = 1\n");
    Unigrid = 1;
  }
*/
 
  /* Initialize the under subgrid field for particle movement. */
 
  if (debug) printf("Call ZeroSUS on TopGrid\n");
 
  Grid->GridData->ZeroSolutionUnderSubgrid(NULL, ZERO_UNDER_SUBGRID_FIELD);
 
  /* Generate this many grids (on this processor). */
 
/*
  Unigrid = 0;
  if (debug) printf("Re-set Unigrid = 0\n");
*/
 
  if (debug) printf("Grid structure: %"ISYM"\n", (int) (sizeof(grid)));
  if (debug) printf("SubGrids structure: %"ISYM"\n", (int) ((Layout[0]*Layout[1]*Layout[2])*sizeof(grid)));
 
  grid *NewGrid, *OldGrid = Grid->GridData;
  grid **SubGrids = new grid*[Layout[0]*Layout[1]*Layout[2]];
  HierarchyEntry *ThisGrid;
 
  int gridcounter = 0;
 
  for (k = 0; k < Layout[2]; k++)
    for (j = 0; j < Layout[1]; j++)
      for (i = 0; i < Layout[0]; i++) {
 
	/* Allocate a new grid hierarchy entry and insert into linked list. */
 
	if (gridcounter == 0)
	  ThisGrid = Grid;
	else {
	  ThisGrid = new HierarchyEntry;
	  ThisGrid->NextGridThisLevel = Grid->NextGridThisLevel;
	  ThisGrid->NextGridNextLevel = NULL;
	  ThisGrid->ParentGrid        = Grid->ParentGrid;
	  Grid->NextGridThisLevel     = ThisGrid;
	}
 
	/* Allocate a new grid and prepare it. */
 
	NewGrid = new grid;
	ThisGrid->GridData = NewGrid;
	NewGrid->InheritProperties(OldGrid);
	NewGrid->SetGravityParameters(OldGrid->ReturnGravityBoundaryType());
 
	/* Compute grid region. */
 
//      printf("GC K J I: %"ISYM" %"ISYM" %"ISYM" %"ISYM"\n",gridcounter,k,j,i);
 
	for (dim = 0; dim < MAX_DIMENSION; dim++) {
	  ijk = (dim == 0) ? i : ((dim == 1) ? j : k);
	  TempDims[dim] = GridDims[dim][ijk];
	  LeftEdge[dim] = Left[dim] + (Right[dim] - Left[dim])*
	    FLOAT(StartIndex[dim][ijk])/FLOAT(Dims[dim]);
	  RightEdge[dim] = Left[dim] + (Right[dim] - Left[dim])*
	    FLOAT(StartIndex[dim][ijk]+TempDims[dim])/FLOAT(Dims[dim]);
	  if (dim < Rank)
	    TempDims[dim] += 2*DEFAULT_GHOST_ZONES;
 
//        printf("  LeftEdge[%"ISYM"] = %8.4"FSYM"  RightEdge[%"ISYM"] = %8.4"FSYM"\n",
//               dim, LeftEdge[dim], dim, RightEdge[dim]);
 
	}
 
	NewGrid->PrepareGrid(Rank, TempDims, LeftEdge, RightEdge, 0);
 
	/* Record this subgrid number in the oldgrid's undersubgrid field. */
 
//      printf("Call ZeroSUS on OldGrid with Value = %10.4e\n", float(gridcounter+1));
 
	if (OldGrid->ZeroSolutionUnderSubgrid(NewGrid,
		   ZERO_UNDER_SUBGRID_FIELD, float(gridcounter+1)) == FAIL) {
	  fprintf(stderr, "Error in grid->ZeroSolutionUnderSubgrid.\n");
	  return FAIL;
	}
	SubGrids[gridcounter] = NewGrid;
 
	gridcounter++;
 
      }
 
 
  Unigrid = 0;
  if (debug) printf("Re-set Unigrid = 0\n");
 
  /* Move Particles (while still on same processor). */
 
  if (OldGrid->MoveSubgridParticlesFast(gridcounter, SubGrids, TRUE) == FAIL) {
    fprintf(stderr, "Error in grid->MoveSubgridParticlesFast.\n");
    return FAIL;
  }
 
  delete [] SubGrids;
 
  /* Distribute new grids amoung processors (and copy out fields). */

#ifdef USE_MPI 
  MPI_Barrier(MPI_COMM_WORLD);
#endif /* USE_MPI */
 
  gridcounter = 0;
  ThisGrid = Grid;

  if (MyProcessorNumber == ROOT_PROCESSOR) 
    printf("Grid distribution\n");
 
  for (k = 0; k < Layout[2]; k++)
    for (j = 0; j < Layout[1]; j++)
      for (i = 0; i < Layout[0]; i++) {
 
	grid *NewGrid = ThisGrid->GridData;
 
	/* Broadcast the number of particles to the other processors
	   (OldGrid is assumed to be on the root processor). */
 
	if (NumberOfProcessors > 1) {

	  int IntTemp = NewGrid->ReturnNumberOfParticles();
 
//          printf("NewGrid->ReturnNumberOfParticles: %"ISYM"\n", IntTemp);
 
	  CommunicationBroadcastValue(&IntTemp, ROOT_PROCESSOR);

	  NewGrid->SetNumberOfParticles(IntTemp);

//          printf("NG particle number set to %"ISYM"\n", IntTemp);

	}
 
	/* Transfer from Old to New (which is still also on root processor) */
 
	FLOAT Zero[] = {0,0,0};
 
	if (MyProcessorNumber == ROOT_PROCESSOR)
	  NewGrid->AllocateGrids();
 
        if (PartitionNestedGrids == FALSE) {
          if (NewGrid->CopyZonesFromGrid(OldGrid, Zero) == FAIL) {
            fprintf(stderr, "Error in grid->CopyZonesFromGrid.\n");
            return FAIL;
          }
        }
 
 
	/* Set processor number of new grid.  Cyclic distribution. */
 
        int NewProc = gridcounter % NumberOfProcessors;
        int ProcMap = ABS(NewProc - NumberOfProcessors) % NumberOfProcessors;
 
        if(NewGrid->ReturnGridInfo(&Rank, Dims, LeftEdge, RightEdge) == FAIL) {
          fprintf(stderr, "Error in grid->ReturnGridInfo.\n");
          return FAIL;
        }
 
	/* Move Grid from current processor to new Processor. */
 
#ifdef USE_OLD_CPU_DISTRIBUTION
	NewGrid->CommunicationMoveGrid(NewProc);
        if (MyProcessorNumber == ROOT_PROCESSOR) {
          printf("Grid = %"ISYM", K J I: [%"ISYM",%"ISYM",%"ISYM"] Proc = %"ISYM"\n", gridcounter, k, j, i, NewProc);
          for (dim = 0; dim < Rank; dim++) {
            printf("  %"ISYM" ::  LeftEdge[%"ISYM"] = %8.4"PSYM"  RightEdge[%"ISYM"] = %8.4"PSYM"\n",
                   NewProc, dim, LeftEdge[dim], dim, RightEdge[dim]);
          }
        }

#endif

#ifdef USE_PERMUTED_CPU_DISTRIBUTION
	NewGrid->CommunicationMoveGrid(ProcMap);
        if (MyProcessorNumber == ROOT_PROCESSOR) {
          printf("Grid = %"ISYM", K J I: [%"ISYM",%"ISYM",%"ISYM"] Proc = %"ISYM"\n", gridcounter, k, j, i, ProcMap);
          for (dim = 0; dim < Rank; dim++) {
            printf("  %"ISYM" ::  LeftEdge[%"ISYM"] = %8.4"PSYM"  RightEdge[%"ISYM"] = %8.4"PSYM"\n",
                   ProcMap, dim, LeftEdge[dim], dim, RightEdge[dim]);
          }
        }
#endif

        /* Detach ForcingFields from BaryonFields. */
 
        if (RandomForcing == 1 && ParallelRootGridIO != 1)
          NewGrid->RemoveForcingFromBaryonFields(); //AK
 
	gridcounter++;
	ThisGrid = ThisGrid->NextGridThisLevel;
 
      }

#ifdef USE_MPI
  MPI_Barrier(MPI_COMM_WORLD);
#endif /* USE_MPI */
 
  /* Clean up. */

  if (RandomForcing == 1 && ParallelRootGridIO != 1)
    OldGrid->RemoveForcingFromBaryonFields();
 
  if (MyProcessorNumber == ROOT_PROCESSOR)
    printf("Delete OldGrid\n");
 
  delete OldGrid;
 
  if (MyProcessorNumber == ROOT_PROCESSOR)
    printf("OldGrid deleted\n");
 
  for (dim = 0; dim < MAX_DIMENSION; dim++) {
    delete [] GridDims[dim];
    delete [] StartIndex[dim];
  }
 
  //  printf("Exit CommunicationPartitionGrid on CPU %"ISYM"\n", MyProcessorNumber);
#ifdef USE_MPI
  MPI_Barrier(MPI_COMM_WORLD);
#endif /* USE_MPI */
 
  return SUCCESS;
}
