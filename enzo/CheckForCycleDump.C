/***********************************************************************
/
/  ENZO Test for CYCLE file in cwd
/
/  written by: Robert Harkness
/  date:       April, 2010
/
************************************************************************/

#include <stdio.h>

#ifdef USE_MPI
#include <mpi.h>
#endif

#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "Hierarchy.h"
#include "TopGridData.h"
#include "LevelHierarchy.h"
#include "CosmologyParameters.h"
#include "version.def"


int CommunicationBroadcastValue(int *Value, int BroadcastProcessor);


int CheckForCycleDump(void)
{

  FILE *con;
  int flag, zero;

  flag = 1000000;
  zero = 0;

  if (MyProcessorNumber == 0) {

    con = fopen("CYCLE", "r");

    if (con != NULL) {
      fscanf(con, "%"ISYM, &flag);
      fclose(con);
    }

    fprintf(stderr, "CycleFlag = %"ISYM"\n", flag);

  }

#ifdef USE_MPI
  CommunicationBroadcastValue(&flag, zero);
#endif

  return (flag);

}
