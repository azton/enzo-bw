/***********************************************************************
/
/  ENZO Test for DUMP file in cwd
/
/  written by: Robert Harkness
/  date:       January, 2010
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
void my_exit(int status);

int CheckForDump(void)
{

  FILE *con;
  int flag, zero;

  flag = 0;
  zero = 0;

  if (MyProcessorNumber == 0) {

    con = fopen("DUMP", "r");

    if (con != NULL) {
      flag = 1;
      fclose(con);
    }

    fprintf(stderr, "DumpFlag = %"ISYM"\n", flag);

  }

  CommunicationBroadcastValue(&flag, zero);
  return (flag);

}
