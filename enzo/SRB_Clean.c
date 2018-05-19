#include <stdio.h>
#include <string.h>
#include <stdlib.h>

int SRB_Delete(char* srb_directory, char* srb_file);




int SRB_Clean(char* filename, char* srbdir)
{
  int status;
  char *srbname;
  char *i;

  fprintf(stderr, "SRB_Clean: %s\n", filename);

  srbname = malloc(32);

  i = rindex(filename, '/');
  if ( i != NULL )
    strcpy(srbname, i+1);

//  strcpy(srbname, filename);
 
  fprintf(stderr, "SRB_Clean: %s \n", srbname);

//  status = SRB_Delete("/home/harkness.teragrid/NCSA", srbname);

  status = SRB_Delete(srbdir, srbname);
  fprintf(stderr,"Delete status: %d\n", status);

  if ( status < 0 )
    status = 0;

  return (status);
}
