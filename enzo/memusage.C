#include <stdio.h>
#include <string.h>

long long int memusage(void)
{
  long long int mem;
  long long kb = 1024;

  FILE *readmem;
  char line[128];
  readmem = fopen("/proc/self/status", "r");

  while(fgets(line,128,readmem) != NULL)
  {
    if (strncmp(line, "VmSize:", 7) == 0)
      sscanf(line, "VmSize: %lld", &mem);
  }
  fclose(readmem);

  mem = mem * kb;
  return mem;
}

