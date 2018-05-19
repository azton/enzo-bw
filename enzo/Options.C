/***********************************************************************
/
/  ENZO Version and Options in Effect
/
/  written by: Robert Harkness
/  date:       15th December 2010
/
************************************************************************/
 
#include <stdio.h>
#include <strings.h>
 
#ifdef USE_MPI
#include <mpi.h>
#endif

#ifdef HYBRID
#include <omp.h>
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
#include "fortran.def"
 
#ifdef MEM_TRACE
Eint64 mused(void);
#endif
 
 
int ENZO_OptionsinEffect(void) 
{

  FILE *opf;

  if (MyProcessorNumber == 0) {

    opf = fopen("Enzo_Options", "w");

    fprintf(opf, "ENZO Options in Effect\n");

/*
    fprintf(opf,"ENZO CVS VERSION %d.%d.%d -- %s\n\n",
	   ENZO_VERSION_MAJOR,
	   ENZO_VERSION_MINOR,
	   ENZO_VERSION_REVISION,
           ENZO_VERSION_DATE);
*/

#ifdef SMALL_INTS
    fprintf(opf, " 32 bit Integer version\n");
#endif

#ifdef LARGE_INTS
    fprintf(opf, " 64 bit Integer version\n");
#endif

#ifdef INITS32
    fprintf(opf, " 32 bit Integer initial conditions\n");
#endif

#ifdef INITS64
    fprintf(opf, " 64 bit Integer initial conditions\n");
#endif

#ifdef r4
    fprintf(opf, " Float precision is 32 bits\n");
#endif

#ifdef r8
    fprintf(opf, " Float precision is 64 bits\n");
#endif

#ifdef p4
    fprintf(opf, " Position and time precision is 32 bits - NOT SUPPORTED!\n");
#endif

#ifdef p8
    fprintf(opf, " Position and time precision is 64 bits\n");
#endif

#ifdef p16
    fprintf(opf, " Position and time precision is 128 bits\n");
#endif

#ifdef HDF5_BE
    fprintf(opf, " HDF5 native i/o type is big-endian\n");
#endif

#ifdef HDF5_LE
    fprintd(opf, " HDF5 native i/o type is little-endian\n");
#endif

#ifdef USE_MPI
    fprintf(opf, " MPI    NumberOfTasks: %d\n", (Eint32) NumberOfProcessors);
#endif

#ifdef HYBRID

    Eint32 NumberOfThreads;
    Eint32 NumberOfProcs;

#pragma omp parallel private(NumberOfThreads) shared(opf) default(none)
  {
    NumberOfThreads = omp_get_max_threads();
#pragma omp master
    fprintf(opf, " OpenMP NumberOfThreads: %d\n", (Eint32) NumberOfThreads);
  }

    NumberOfProcs = omp_get_num_procs();
    fprintf(opf, " OpenMP NumberOfCores: %d\n", (Eint32) NumberOfProcs);

#else

    fprintf(opf, " No OpenMP Threading\n");

#endif


    fprintf(opf, "\n");
    fprintf(opf, "Optimizations in Effect\n");

#ifdef OOC_BOUNDARY
    fprintf(opf, "  Out-of-core Top Grid boundary conditions\n");
#endif

#ifdef SIB1
    fprintf(opf, "  Fast Sibling Locator 1\n");
#endif

#ifdef SIB2
    fprintf(opf, "  Fast Sibling Locator 2\n");
#endif

#ifdef SIB3
    fprintf(opf, "  Fast Sibling Locator 3\n");
#endif

#ifdef SIB4
    fprintf(opf, "  Fast Sibling Locator 4\n");
#endif

#ifdef SIB5
    fprintf(opf, "  Fast Sibling Locator 5\n");
#endif

#ifdef STATIC_SIBLING_LIST
    fprintf(opf, "  Static allocation of Level Zero Sibling List\n");
#endif

#ifdef FLUX_FIX
    fprintf(opf, "  New Flux Correction scheme by Collins & Wagner\n");
#endif

#ifdef SAB
    fprintf(opf, "  AccelerationHack by Collins\n");
#endif

#ifdef USE_DT_LIMIT
    fprintf(opf, "  Use dt limit in AMR\n");
#endif

#ifdef UNIGRID
    fprintf(opf, "  Minimum memory start-up => non-apative mesh only\n");
#endif

#ifdef UNIGRID_TRANSPOSE
    fprintf(opf, "  Fast book-keeping for FFT TopGrid neighbours\n");
#endif

#ifdef HDF5_USE_HDF5_GROUPS
    fprintf(opf, "  HDF5 groups for packed AMR\n");
#endif

#ifdef USE_HDF5_OUTPUT_BUFFERING
    fprintf(opf, "  HDF5 in-core buffering for packed AMR output\n");
#endif

#ifdef USE_HDF5_INPUT_BUFFERING
    fprintf(opf, "  HDF5 in-core buffering for packed AMR input\n");
#endif

#ifdef SINGLE_HDF5_OPEN_ON_INPUT
    fprintf(opf, "  HDF5 single open on input - no explicit task map\n");
#endif

#ifdef USE_NEW_RESTART_FORMAT
    fprintf(opf, "  New hierarchy format\n")
#endif

#ifdef FORCE_BUFFER_PURGE
    fprintf(opf, "  Force purge of communication buffers\n");
#endif

#ifdef FORCE_MSG_PROGRESS
    fprintf(opf, "  Force message progress with MPI_Barrier calls\n");
#endif

#ifdef ENABLE_LOAD_BALANCE
    fprintf(opf, "  Load balancing enabled\n");
#else
    fprintf(opf, "  Load balancing disabled\n");
#endif

#ifdef ISOLATED_GRAVITY
    fprintf(opf, "  Isolated gravity enabled\n");
#else /* ISOLATED_GRAVITY */
    fprintf(opf, "  Isolated gravity disabled\n");
#endif /* ISOLATED_GRAVITY */

#ifdef RATE_AND_COOL
    fprintf(opf, "  Dual solve rate and cool enabled\n");
#endif

#ifdef MEM_TRACE
    fprintf(opf, "  Memory tracing enabled\n");
#endif

#ifdef TP_VELOCITY
    fprintf(opf, "  Tracer particle velocity output is enabled\n");
#else
    fprintf(opf, "  Tracer particle velocity output is disabled\n");
#endif

#ifdef RAD_HYDRO
    fprintf(opf, "  RHD is enabled\n");
#else
    fprintf(opf, "  RHD is disabled\n");
#endif

#ifdef USE_CRAYX1_FFT
    fprintf(opf, "Cray X1 fft\n");
#endif
#ifdef USE_ACML_FFT
    fprintf(opf, "ACML fft\n");
#endif
#ifdef USE_MKL_FFT
    fprintf(opf, "MKL fft\n");
#endif
#ifdef USE_ESSL_FFT
    fprintf(opf, "ESSL fft\n");
#endif
#ifdef USE_SGI_FFT
    fprintf(opf, "SGI SCSL fft\n");
#endif
#ifdef USE_FFTE_FFT
    fprintf(opf, "FFTE fft\n");
#endif
#ifdef USE_S90_FFT
    fprintf(opf, "Singleton F90 fft\n");
#endif
#ifdef USE_S66_FFT
    fprintf(opf, "Singleton F66 fft\n");
#endif

    fprintf(opf, "\n");
    fprintf(opf, "Macro and Parameter Definitions\n");

    fprintf(opf, "  MAX_NUMBER_OF_TASKS                 %8d\n", MAX_NUMBER_OF_TASKS);
    fprintf(opf, "  MAX_NUMBER_OF_NODES                 %8d\n", MAX_NUMBER_OF_NODES);
    fprintf(opf, "  MAX_TASKS_PER_NODE                  %8d\n", MAX_TASKS_PER_NODE);
    fprintf(opf, "  MAX_NUMBER_OF_FIELD_TYPES,          %8d\n", MAX_NUMBER_OF_FIELD_TYPES);
    fprintf(opf, "  MAX_NUMBER_OF_BARYON_FIELDS (>=6)   %8d\n", MAX_NUMBER_OF_BARYON_FIELDS);
    fprintf(opf, "  MAX_NUMBER_OF_SUBGRIDS              %8d\n", MAX_NUMBER_OF_SUBGRIDS);
    fprintf(opf, "  MAX_DEPTH_OF_HIERARCHY              %8d\n", MAX_DEPTH_OF_HIERARCHY);
    fprintf(opf, "  MAX_LINE_LENGTH                     %8d\n", MAX_LINE_LENGTH);

    fprintf(opf, "  MAX_NAME_LENGTH                     %8d\n", MAX_NAME_LENGTH);
    fprintf(opf, "  MAX_GRID_TAG_SIZE                   %8d\n", MAX_GRID_TAG_SIZE);
    fprintf(opf, "  MAX_TASK_TAG_SIZE                   %8d\n", MAX_TASK_TAG_SIZE);
    fprintf(opf, "  MAX_GROUP_TAG_SIZE                  %8d\n", MAX_GROUP_TAG_SIZE);
    fprintf(opf, "  MAX_CYCLE_TAG_SIZE                  %8d\n", MAX_CYCLE_TAG_SIZE);

    fprintf(opf, "  GRID FORMAT                              "); fprintf(opf, GRID_TAG_FORMAT);  fprintf(opf, "\n");
    fprintf(opf, "  TASK FORMAT                              "); fprintf(opf, TASK_TAG_FORMAT);  fprintf(opf, "\n");
    fprintf(opf, "  GROUP FORMAT                             "); fprintf(opf, GROUP_TAG_FORMAT); fprintf(opf, "\n");
    fprintf(opf, "  CYCLE FORMAT                             "); fprintf(opf, CYCLE_TAG_FORMAT); fprintf(opf, "\n");

    fprintf(opf, "  MAX_COUNTERS                        %8d\n", MAX_COUNTERS);
    fprintf(opf, "  DEFAULT_GHOST_ZONES (>=3)           %8d\n", DEFAULT_GHOST_ZONES);
    fprintf(opf, "  MAX_NUMBER_OF_OUTPUT_REDSHIFTS      %8d\n", MAX_NUMBER_OF_OUTPUT_REDSHIFTS);
    fprintf(opf, "  GRAVITY_BUFFER_SIZE                 %8d\n", GRAVITY_BUFFER_SIZE);
    fprintf(opf, "  MAX_FLAGGING_METHODS                %8d\n", MAX_FLAGGING_METHODS);
    fprintf(opf, "  MAX_STATIC_REGIONS                  %8d\n", MAX_STATIC_REGIONS);
    fprintf(opf, "  MAX_NUMBER_OF_PARTICLE_ATTRIBUTES   %8d\n", MAX_NUMBER_OF_PARTICLE_ATTRIBUTES);
    fprintf(opf, "  MAX_TIME_ACTIONS                    %8d\n", MAX_TIME_ACTIONS);
    fprintf(opf, "  MAX_CUBE_DUMPS                      %8d\n", MAX_CUBE_DUMPS);
    fprintf(opf, "  MAX_POTENTIAL_ITERATIONS            %8d\n", MAX_POTENTIAL_ITERATIONS);
    fprintf(opf, "  MAX_COLOR                           %8d\n", MAX_COLOR);
    fprintf(opf, "  MAX_ANY_SINGLE_DIRECTION            %8d\n", MAX_ANY_SINGLE_DIRECTION);

    fprintf(opf, "\n");

#ifdef LINUX_PGI_X86_64
    fprintf(opf, "Linux PGI x86_64\n");
#endif

#ifdef LINUX_INTEL_X86_64
    fprintf(opf, "Linux Intel x86_64\n");
#endif

#ifdef LINUX_IBM_XL
    fprintf(opf, "Linux IBM XL\n");
#endif

#ifdef AIX_IBM_XL
    fprintf(opf, "AIX IBM XL\n");
#endif

#ifdef LINUX_CRAY_X86_64
    fprintf(opf, "Linux Cray CCE x86_64\n");
#endif

#ifdef LINUX_CRAY_X2
    fprintf(opf, "Linux Cray X2\n");
#endif

#ifdef XT3
    fprintf(opf, "Processor type is CRAY XT\n");
#endif

    fclose(opf);

  } // processor zero only

  return SUCCESS;

}
