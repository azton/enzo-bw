#include <stdio.h>
#include <string.h>
#include <stdlib.h>

int SRB_Put(char *srb_directory, char *srb_file, char *local_file);

int SRB_Mover(char* filename, char* srbdir)
{
  int status;
  char *srbname;
  char *i;

  fprintf(stderr, "SRB_Put: %s\n", filename);

  srbname = malloc(32);

  i = rindex(filename, '/');
  if ( i != NULL )
    strcpy(srbname, i+1);

//  strcpy(srbname, filename);
 
  fprintf(stderr, "SRB_File: %s ==> %s/%s\n", filename, srbdir, srbname);

//  status = srb_file_check("/home/harkness.teragrid/NCSA", srbname);
//  fprintf(stderr,"Delete status: %d\n", status);

//  status = SRB_Put("/home/harkness.teragrid/NCSA", srbname, filename);

  status = SRB_Put(srbdir, srbname, filename);
  fprintf(stderr,"Put status: %d\n", status);

  status = 0;
  return (status);
}
