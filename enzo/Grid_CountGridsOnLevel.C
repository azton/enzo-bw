#include <stdio.h>

#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"

int grid::CountGridsOnLevel(int level)
{

  if (ProcessorNumber != MyProcessorNumber)
    return SUCCESS;

  level_count[level]++;

  return SUCCESS;

} 
