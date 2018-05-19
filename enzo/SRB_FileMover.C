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
 
extern "C" int SRB_Mover(char* file, char* srbdir);
extern "C" int SRB_Clean(char* file, char* srbdir);
 
int FileMover( char *file )
{
  int status;
  char srbdir[120];
 
//  strcpy(srbdir, SRBprefix);
 
  strcpy(srbdir, SRBcwd);
 
  status = SRB_Clean(file, srbdir);
  fprintf(stderr, "FileMover - Clean: %s  status %"ISYM"\n", file, status);
    assert( status == 0 );
 
  status = SRB_Mover(file, srbdir);
  fprintf(stderr, "FileMover - Mover: %s  status %"ISYM"\n", file, status);
    assert( status == 0 );
 
  return (status);
}
