#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
 
#ifdef USE_MPI
#include "mpi.h"
#endif /* USE_MPI */
 
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
 
 
int MakeDir( char *subdirectory )
{
  char srbdir[120];
 
  strcpy(srbdir, SRBprefix);
 
  fprintf(stderr, "Dummy MakeDir %"ISYM": %s / %s\n", MyProcessorNumber, srbdir, subdirectory);
 
  strcpy(SRBcwd, SRBprefix);
  strcat(SRBcwd, "/");
  strcat(SRBcwd, subdirectory);
 
  fprintf(stderr, "SRB cwd %"ISYM": %s\n", MyProcessorNumber, SRBcwd);
 
  return (0);
}
