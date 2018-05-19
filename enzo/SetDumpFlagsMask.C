/***********************************************************************
/
/  SET DATA DUMP FLAGS MASK
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


int SetDumpFlagsMask(void)
{

  int i;

  DumpTemperature = 0;
  DumpDarkMatterDensity = 0;
  DumpParticles = 0;

  for ( i = 0; i < MAX_NUMBER_OF_FIELD_TYPES; i++ ) {
    DumpBaryonField[i] = 0;
  }

  // Should read in a set of flags for this

  /*
  Density         = 0,
  TotalEnergy     = 1,
  InternalEnergy  = 2,
  Pressure        = 3,
  Velocity1       = 4,
  Velocity2       = 5,
  Velocity3       = 6,
  ElectronDensity = 7,
  HIDensity       = 8,
  HIIDensity      = 9,
  HeIDensity      = 10,
  HeIIDensity     = 11,
  HeIIIDensity    = 12,
  */

  DumpBaryonField[Density] = 1;

  return SUCCESS;

}
