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
 
extern "C" int SRB_Mkdir(char *srb_directory, char *srb_subdir);
 
int MakeDir( char *subdirectory )
{
  int status;
  char srbdir[40];
 
  strcpy(srbdir, SRBprefix);
 
  status = SRB_Mkdir(srbdir, subdirectory);
  fprintf(stderr, "MakeDir: %s %s status %"ISYM"\n", srbdir, subdirectory, status);
    assert( status == 0 );
 
  strcpy(SRBcwd, SRBprefix);
  strcat(SRBcwd, "/");
  strcat(SRBcwd, subdirectory);
 
  fprintf(stderr, "SRB cwd: %s\n", SRBcwd);
 
  return (status);
}
