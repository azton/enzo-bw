/***********************************************************************
/
/  GRID CLASS (FLAG CELLS TO BE REFINED BY MASS)
/
/  written by: Greg Bryan
/  date:       November, 1994
/  modified1:
/
/  PURPOSE:
/
/  RETURNS:
/    number of flagged cells, or -1 on failure
/
************************************************************************/
 
#include <stdio.h>
#include <math.h>
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
 
int grid::FlagCellsToBeRefinedByMass(int level, int method)
{
 
  /* Return if this grid is not on this processor. */
 
  if (MyProcessorNumber != ProcessorNumber)
    return SUCCESS;
 
  /* error check */
 
  if (MassFlaggingField == NULL || FlaggingField == NULL) {
    fprintf(stderr, "MassFlaggingField or Flagging Field is undefined.\n");
    return -1;
  }
 
  /* compute size */
 
  int i, size = 1;
  for (int dim = 0; dim < GridRank; dim++)
    size *= GridDimension[dim];
 
  /* Compute the ModifiedMinimumMass */
 
  float ModifiedMinimumMassForRefinement =
    MinimumMassForRefinement[method]*POW(RefineBy,
		    level*MinimumMassForRefinementLevelExponent[method]);
  if (ProblemType == 28)
    ModifiedMinimumMassForRefinement = 0;
 
  /* Flag points */
 
  for (i = 0; i < size; i++)
    FlaggingField[i] += (MassFlaggingField[i]
			     > ModifiedMinimumMassForRefinement ? 1 : 0);
 
  /* Count number of flagged Cells. */
 
  int NumberOfFlaggedCells = 0;
  for (i = 0; i < size; i++) {
    FlaggingField[i] = min(FlaggingField[i], 1);
    NumberOfFlaggedCells += FlaggingField[i];
  }
 
  /* remove MassFlaggingField. */
 
  delete MassFlaggingField;
  MassFlaggingField = NULL;
 
  return NumberOfFlaggedCells;
 
}
