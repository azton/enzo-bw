/***********************************************************************
/
/  INITIALIZE LOCAL COMPONENTS OF A NEW SIMULATION
/
/  written by: Daniel Reynolds
/  date:       April 2006
/
/  PURPOSE:
/
/  RETURNS: SUCCESS or FAIL
/
************************************************************************/
 
// This routine intializes the local components of a new simulation 

#ifdef RAD_HYDRO
#include "ImplicitProblemABC_preincludes.h"
#endif
 
#include <string.h>
#include <stdio.h>
 
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "Hierarchy.h"
#include "TopGridData.h"
#include "StarParticleData.h"
#ifdef RAD_HYDRO
#include "ImplicitProblemABC.h"
#endif
 
// Function prototypes
 
int CosmologySimulationReInitialize(HierarchyEntry *TopGrid,
                                    TopGridData &MetaData);
 
int NestedCosmologySimulationReInitialize(HierarchyEntry *TopGrid,
                                          TopGridData &MetaData);
 
int TurbulenceSimulationReInitialize(HierarchyEntry *TopGrid,
                                    TopGridData &MetaData);
 
 
 

int InitializeLocal(int restart, HierarchyEntry &TopGrid, TopGridData &MetaData
#ifdef RAD_HYDRO
		    , ImplicitProblemABC *ImplicitSolver
#endif
)
{

  // If performing a radiation-hydrodynamics simulation (200 series),
  // set up the ImplicitSolver problem based on grid structures, etc.
  // (assumes that DetermineParallelism has already been run on the top grid)
  //   if ((ProblemType > 199) && (ProblemType < 300)) {

#ifdef RAD_HYDRO
  if (RadiationHydrodynamics > 0) {
    if (ImplicitSolver->Initialize(TopGrid, MetaData) == FAIL) {
      fprintf(stderr,"Error Initializing ImplicitSolver solver\n");
      return FAIL;
    }
  }
#endif

  // Call local problem initializer
  if (debug)
    printf("InitializeLocal: Starting problem initialization.\n");

  // For problem 30 if starting from scratch, using ParallelGridIO,
  // read in data only after partitioning the grid

/* 
  if (!restart) {
    if (debug)
      if (ParallelRootGridIO == TRUE && ProblemType == 30) {
	if (PartitionNestedGrids) {
	  printf("InitializeLocal: Re-initialize NestedCosmologySimulation\n");
	} else {
	  printf("InitializeLocal: Re-initialize CosmologySimulation\n");
	}
      }
    
    if (ParallelRootGridIO == TRUE && ProblemType == 30) {
      if (PartitionNestedGrids) {
	if (NestedCosmologySimulationReInitialize(&TopGrid, MetaData) == FAIL) {
	  fprintf(stderr, "Error in NestedCosmologySimulationReInitialize.\n");
	  return FAIL;
	}
      } else {
	if (CosmologySimulationReInitialize(&TopGrid, MetaData) == FAIL) {
	  fprintf(stderr, "Error in CosmologySimulationReInitialize.\n");
	  return FAIL;
	}
      }
    }
  }
 
  // For problem 60 if starting from scratch, using ParallelGridIO, 
  // read in data only after partitioning grid.

  if (!restart) {
    if (ParallelRootGridIO == TRUE && ProblemType == 60)
      if (TurbulenceSimulationReInitialize(&TopGrid, MetaData) == FAIL) {
	fprintf(stderr, "Error in TurbulenceSimulationReInitialize.\n");
	return FAIL;
      }
  }
*/

  // Insert new problem intializer here...
 
  if (debug)
    printf("InitializeLocal: Finished problem initialization.\n");
 

  return SUCCESS;
 
}
