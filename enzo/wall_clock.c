#include <stdio.h>
#include <time.h>
#include <sys/times.h>
#include <sys/time.h>

#if defined(LINUX_PGI_X86_64) || defined(LINUX_INTEL_X86_64) || defined(LINUX_CRAY_X86_64) || defined(LINUX_CRAY_X2) || defined(XT3)
double wall_clock_(void);

double wall_clock_(void)
#endif
#if defined(LINUX_IBM_XL) || defined(AIX_IBM_XL)
double wall_clock(void);

double wall_clock(void)
#endif
{
  struct timeval timestr;
  void *Tzp=0;
  gettimeofday(&timestr, Tzp);
  return (double)timestr.tv_sec+1.0E-06*(double)timestr.tv_usec;
}
