/***********************************************************************
/
/  UPDATE SOFT LINK TO LAST RESTART DUMP
/
/  written by: Robert Harkness
/  date:       May, 2008
/
/  PURPOSE:
/
************************************************************************/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <assert.h>

#ifdef USE_MPI
#include "mpi.h"
#endif /* USE_MPI */

#include "macros_and_parameters.h"


int UpdateSymbolicLink(char *GlobalDir)
{

  FILE *fptr;

//  char GlobalDir[MAX_LINE_LENGTH];
  char DummyLink[MAX_LINE_LENGTH];
  char line[MAX_LINE_LENGTH];
  char LastDump[MAX_LINE_LENGTH];
  int status;

  LastDump[0] = 0;

//  Find the last SRBlog entry

  fptr = fopen("SRBlog", "r");
  while (fgets(line, MAX_LINE_LENGTH, fptr) != NULL) {
    sscanf(line, "%s", LastDump);
  }
  fclose(fptr);

  fprintf(stderr, "%s\n", LastDump);

//  Link this file to a dummy filename in GlobalDir

//  strcpy(GlobalDir, "/scratch/00770/harkness/LAF/Dumps");
//  fprintf(stderr, "%s\n", GlobalDir);

  strcpy(DummyLink, GlobalDir);
  strcat(DummyLink, "/ENZO_Restart");
  fprintf(stderr, "%s\n", DummyLink);

//  Unlink the dummy file!

  status = 0;
  status = unlink(DummyLink);
  fprintf(stderr, "Unlink status %d\n", status);


  status = 0;
  status = symlink(LastDump, DummyLink);
  fprintf(stderr, "Symlink status %d\n", status);

  strcat(DummyLink, ".hierarchy");
  strcat(LastDump, ".hierarchy");

  status = 0;
  status = unlink(DummyLink);
  fprintf(stderr, "Unlink status %d\n", status);


  status = 0;
  status = symlink(LastDump, DummyLink);
  fprintf(stderr, "Symlink status %d\n", status);

  return 0;
}
