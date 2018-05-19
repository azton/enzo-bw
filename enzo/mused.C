#if defined(AIX_IBM_XL)

// Memory usage from getrusage

#include<stdio.h>
#include<sys/time.h>
#include<sys/resource.h>
#include<unistd.h>
#include<assert.h>

long long int mused(void)
{

#ifdef MEM_TRACE
  struct rusage temp;
  long long int bytes;
  int result;

  result = getrusage(RUSAGE_SELF, &temp);
  if( result == 0 ) {
    bytes = ((long long int) (1024)) * ((long long int) temp.ru_maxrss);
  } else {
    bytes = ((long long int) (0));
  }
  return(bytes);
#else
  return((long long int) 0);
#endif

}

#elif defined(LINUX_PGI_X86_64) || defined(LINUX_INTEL_X86_64) || defined(LINUX_IBM_XL) || defined(LINUX_CRAY_X86_64) || defined(XT3)


#include <stdio.h>
#include <string.h>

long long int mused(void)
{

#ifdef MEM_TRACE
  long long int mem;
  long long kb = 1024;

  FILE *readmem;
  char line[128];
  readmem = fopen("/proc/self/status", "r");

  while(fgets(line,128,readmem) != NULL)
  {
    if (strncmp(line, "VmRSS:", 6) == 0)
      sscanf(line, "VmRSS: %lld", &mem);
  }
  fclose(readmem);

  mem = mem * kb;
  return mem;
#else
  return ((long long int) 0);
#endif

}

#else

// Default case return zero bytes

long long int mused(void)
{
  return((long long int) 0);
}

#endif
