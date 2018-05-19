/***********************************************************************
/
/  GRID CLASS (MOVE A GRID FROM ONE PROCESSOR TO ANOTHER)
/
/  written by: Greg Bryan
/  date:       December, 1997
/  modified1:
/
/  PURPOSE:
/
/  INPUTS:
/
************************************************************************/
 
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
 
/* function prototypes */
 
 
 
int grid::CommunicationMoveGrid(int ToProcessor)
{
 
  int Zero[] = {0, 0, 0};
 
  if ((MyProcessorNumber == ProcessorNumber ||
       MyProcessorNumber == ToProcessor) &&
      ProcessorNumber != ToProcessor) {
 
    /* Copy baryons. */
 
    if (NumberOfBaryonFields > 0)
      this->CommunicationSendRegion(this, ToProcessor, ALL_FIELDS,
				    NEW_ONLY, Zero, GridDimension);
 
    /* Copy particles. */
 
    if (NumberOfParticles > 0)
      this->CommunicationSendParticles(this, ToProcessor, 0,
				       NumberOfParticles, 0);
  }
 
  /* Delete fields on old grid. */
 
  if (MyProcessorNumber == ProcessorNumber && ProcessorNumber != ToProcessor)
    this->DeleteAllFields();
 
  /* Update processor number. */
 
  ProcessorNumber = ToProcessor;
 
  return SUCCESS;
}
 
