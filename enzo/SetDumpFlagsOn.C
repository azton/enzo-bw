/***********************************************************************
/
/  SET ALL DATA DUMP FLAGS TO ON
/
/  written by: Robert Harkness
/  date:       June 2009
/
************************************************************************/


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




int SetDumpFlagsOn(void)
{

  int i;

  DumpTemperature = 1;
  DumpDarkMatterDensity = 1;
  DumpParticles = 1;

  for ( i = 0; i < MAX_NUMBER_OF_FIELD_TYPES; i++ ) {
    DumpBaryonField[i] = 1;
  }

  return SUCCESS;

}
