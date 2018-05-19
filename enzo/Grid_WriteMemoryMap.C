/***********************************************************************
/
/  GRID CLASS (WRITE OUT TASK MEMORY MAP)
/
/  written by: Robert Harkness
/  date:       January 2006
/
/  PURPOSE:
/
************************************************************************/
 
#include <mpi.h>
#include <string.h>
#include <stdio.h>
#include <assert.h>
 
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
 
// function prototypes
 
int grid::WriteMemoryMap(FILE *fptr, char *base_name, int grid_id)
{
 
  int i, j, k;
  int file_status;
  int output_cube;

  if (MyProcessorNumber == ROOT_PROCESSOR) {

    fprintf(fptr, "Grid %8"ISYM"  PN %8"ISYM"  Memory %16"ISYM"\n", grid_id, ProcessorNumber, TaskMemory[ProcessorNumber]);

  }

  return SUCCESS;
 
}
