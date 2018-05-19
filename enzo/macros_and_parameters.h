#ifndef __macros_and_parameters_h_
#define __macros_and_parameters_h_
/***********************************************************************
/  
/ MACRO DEFINITIONS AND PARAMETERS
/
************************************************************************/

#include "message.h"


/* Modifiable Parameters */

#define MAX_NUMBER_OF_TASKS             32768

#define MAX_NUMBER_OF_NODES              1024

#define MAX_TASKS_PER_NODE                   __max_cpu_per_node

#define MAX_NUMBER_OF_BARYON_FIELDS          __max_baryons  /* must be at least 6 */

#define MAX_NUMBER_OF_SUBGRIDS               __max_subgrids

#define MAX_DEPTH_OF_HIERARCHY             40

#define MAX_LINE_LENGTH                   512

#define MAX_NAME_LENGTH                   512

#define MAX_GRID_TAG_SIZE         16
#define MAX_TASK_TAG_SIZE         16
#define MAX_GROUP_TAG_SIZE        16
#define MAX_CYCLE_TAG_SIZE        16

#ifdef LARGE_TAGS
#define GRID_TAG_FORMAT        "8.8"
#define TASK_TAG_FORMAT        "8.8"
#endif
#ifdef SMALL_TAGS
#define GRID_TAG_FORMAT        "4.4"
#define TASK_TAG_FORMAT        "4.4"
#endif

#define GROUP_TAG_FORMAT       "8.8"
#define CYCLE_TAG_FORMAT       "4.4"
#define MAX_COUNTERS              40


#define DEFAULT_GHOST_ZONES                 3  /* at least 3 */

#define MAX_NUMBER_OF_FIELD_TYPES         100

#define MAX_NUMBER_OF_OUTPUT_REDSHIFTS    500

#define GRAVITY_BUFFER_SIZE                 3

#define MAX_FLAGGING_METHODS               13  /* check def in global_data */ 

#define MAX_STATIC_REGIONS               1000

#define MAX_NUMBER_OF_PARTICLE_ATTRIBUTES  10

#define MAX_TIME_ACTIONS                   10

#define MAX_CUBE_DUMPS                     50

#define MAX_CUBE_DUMPS                     50

#define MAX_MOVIE_FIELDS                    1

#define MAX_POTENTIAL_ITERATIONS            0

#define ROOT_PROCESSOR                      0

#define VERSION                           1.3  /* current version number */

/* Unmodifiable Parameters */

#define MAX_DIMENSION                       3  /* must be 3! */

/* Fortran name generator (cpp blues) */

#if defined(SUN_OLD)
#define FORTRAN_NAME(NAME) NAME/**/_
#endif

#if defined(LINUX_PGI_X86_64) || defined(LINUX_INTEL_X86_64) || defined(LINUX_CRAY_X86_64) || defined(LINUX_CRAY_X2) || defined(XT3)
#define FORTRAN_NAME(NAME) NAME##_
#endif

#if defined(LINUX_IBM_XL) || defined(AIX_IBM_XL) || defined(HP_UX)
#define FORTRAN_NAME(NAME) NAME
#endif

/* Precision-related definitions. */

typedef long long long_int;
typedef long double long_double;
typedef unsigned int unsigned_int;
typedef unsigned long long int unsigned_long_int;

/* Previously in hdf4.h */

typedef float        float32;
typedef double       float64;
typedef long double  float128;

/* Macro definitions for portability */

typedef void           *VOIDP;
typedef int            Eint32;
typedef long long int  Eint64;
typedef float          Eflt32;
typedef double         Eflt64;
typedef long double    Eflt128;
typedef long long int  Elong_int;

typedef int            MPI_Arg;

typedef int            HDF5_hid_t;

/* HDF5 definitions */

#ifdef HDF5_BE
#define HDF5_FILE_I4 H5T_STD_I32BE
#define HDF5_FILE_I8 H5T_STD_I64BE
#define HDF5_FILE_R4 H5T_IEEE_F32BE
#define HDF5_FILE_R8 H5T_IEEE_F64BE
#define HDF5_FILE_B8 H5T_STD_B8BE
#endif /* HDF5_BE */

#ifdef HDF5_LE
#define HDF5_FILE_I4 H5T_STD_I32LE
#define HDF5_FILE_I8 H5T_STD_I64LE
#define HDF5_FILE_R4 H5T_IEEE_F32LE
#define HDF5_FILE_R8 H5T_IEEE_F64LE
#define HDF5_FILE_B8 H5T_STD_B8LE
#endif /* HDF5_LE */

#define HDF5_I4  H5T_NATIVE_INT
#define HDF5_I8  H5T_NATIVE_LLONG
#define HDF5_R4  H5T_NATIVE_FLOAT
#define HDF5_R8  H5T_NATIVE_DOUBLE
#define HDF5_R16 H5T_NATIVE_LDOUBLE

/* Precision-dependent definitions */

#if defined(INITS32)
#define inits_type float32
#endif

#if defined(INITS64)
#define inits_type float64
#endif

#define ByteDataType MPI_BYTE

#ifdef SMALL_INTS
#define Eint int
#define Eunsigned_int unsigned_int
#define ISYM "d"
#define LSYM "lld"
#define IntDataType MPI_INT
#define HDF5_INT HDF5_I4
#define HDF5_FILE_INT HDF5_FILE_I4
#define nint(A) ( (int) ((A) + 0.5*sign(A)) )
#define nlongint(A) ( (long_int) ((A) + 0.5*sign(A)) )
#define ABS(A) abs((int) (A))
#endif

#ifdef LARGE_INTS
#define int long_int
#define Eint long_int
#define Eunsigned_int unsigned_long_int
#define ISYM "lld"
#define LSYM "lld"
#define IntDataType MPI_LONG_LONG_INT
#define HDF5_INT HDF5_I8
#define HDF5_FILE_INT HDF5_FILE_I8
#define nint(A) ( (long_int) ((A) + 0.5*sign(A)) )
#define nlongint(A) ( (long_int) ((A) + 0.5*sign(A)) )
#define ABS(A) labs((long_int) (A))
#endif

#ifdef r4
#define Eflt float
#define FSYM "f"
#define ESYM "e"
#define FloatDataType MPI_FLOAT
#ifdef COMPACT_IO
#define HDF5_REAL HDF5_R4
#define HDF5_FILE_REAL HDF5_FILE_R4
#else
#define HDF5_REAL HDF5_R4
#define HDF5_FILE_REAL HDF5_FILE_R8
#endif
#endif

#ifdef r8
#define Eflt double
#define FSYM "lf"
#define ESYM "le"
#define FloatDataType MPI_DOUBLE
#define float32 TEMP_HOLD_NAME
#define float double
#define TEMP_HOLD_NAME float32
#ifdef COMPACT_IO
#define HDF5_REAL HDF5_R8
#define HDF5_FILE_REAL HDF5_FILE_R4
#else
#define HDF5_REAL HDF5_R8
#define HDF5_FILE_REAL HDF5_FILE_R8
#endif
#endif


#ifdef p8
#define FLOAT double
#define PSYM "lf"
#define GSYM "g"
#define GOUTSYM ".14g"
#define MY_MPIFLOAT MPI_DOUBLE
#define FLOATDataType MPI_DOUBLE
#define HDF5_PREC HDF5_R8
#define HDF5_FILE_PREC HDF5_FILE_R8
#endif

#ifdef p16
#define FLOAT long_double
#define PSYM "Lf"
#define GSYM "g"
#define GOUTSYM ".21Lg"
#define MY_MPIFLOAT MPI_LONG_DOUBLE
#define FLOATDataType MPI_LONG_DOUBLE
#define HDF5_PREC HDF5_R16
#define HDF5_FILE_PREC HDF5_FILE_R16
#endif


/* Standard definitions (well, fairly standard) */

#ifndef NULL
#define NULL      0
#endif

#ifdef FAIL
#undef FAIL
#endif
#define FAIL      0
#define SUCCESS   1

#ifndef FALSE
#define FALSE     0
#define TRUE      1
#endif

/* Not-so standard definitions */
#ifndef HDF_FAIL
#define HDF_FAIL -1
#endif

#define FLOAT_UNDEFINED  -99999.0
#define INT_UNDEFINED    -99999

/* Macro definitions (things C should have) */

#define max(A,B) ((A) > (B) ? (A) : (B))
#define min(A,B) ((A) < (B) ? (A) : (B))
#define sign(A)  ((A) >  0  ?  1  : -1 )
#define POW(X,Y) pow((double) (X), (double) (Y))

/* Definitions for FastFourierTransform and related routines */

#define FFT_FORWARD  +1
#define FFT_INVERSE  -1
#define REAL_TO_COMPLEX    0
#define COMPLEX_TO_COMPLEX 1

/* Definitions for grid::RestoreEnergyConsistency */

#define ENTIRE_REGION  0
#define ONLY_BOUNDARY  1

/* Definitions for grid::ZeroSolutionUnderSubgrid */

#define ZERO_ALL_FIELDS          0
#define ZERO_UNDER_SUBGRID_FIELD 1

/* Definitions for grid::CommunicationSend/ReceiveRegion and 
   grid::DepositPositions */

#define MASS_FLAGGING_FIELD              -6
#define ACCELERATION_FIELDS              -5
#define POTENTIAL_FIELD                  -4
#define GRAVITATING_MASS_FIELD           -3
#define GRAVITATING_MASS_FIELD_PARTICLES -2
#define ALL_FIELDS   -1

#define NEW_AND_OLD   0
#define NEW_ONLY      1
#define OLD_ONLY      2

/* Obsolete definitions for grid::ComputeAccelerationField */

#define OBSOLETE_PARTICLES  0
#define OBSOLETE_GRIDS      1
#define OBSOLETE_ZEUS_GRIDS 2

#define DIFFERENCE_TYPE_STAGGERED 0
#define DIFFERENCE_TYPE_NORMAL    1

/* Definitions for CommunicationTranspose */

#define NORMAL_ORDER      0
#define TRANSPOSE_FORWARD 1
#define TRANSPOSE_REVERSE 2

/* Definitions for CommunicationTransferParticles */

#define COPY_IN   0
#define COPY_OUT  1

/* Definitions for CommunicationDirection */

#define COMMUNICATION_SEND_RECEIVE 0
#define COMMUNICATION_RECEIVE      1
#define COMMUNICATION_SEND         2

/* MPI Tags */

#define MPI_TRANSPOSE_TAG 10
#define MPI_SENDREGION_TAG 11
#define MPI_FLUX_TAG 12
#define MPI_TRANSFERPARTICLE_TAG 13
#define MPI_SENDPARTFLOAT_TAG 14
#define MPI_SENDPARTINT_TAG 15

/* Definitions for CommunicationBufferedSend. */

#define BUFFER_IN_PLACE -1

/* Particle types (note: gas is a conceptual type) */

#define NUM_PARTICLE_TYPES 5

#define PARTICLE_TYPE_GAS          0
#define PARTICLE_TYPE_DARK_MATTER  1
#define PARTICLE_TYPE_STAR         2
#define PARTICLE_TYPE_TRACER       3
#define PARTICLE_TYPE_MUST_REFINE  4

#ifdef USE_MPI
#define MPI_INSTRUMENTATION
#endif /* USE_MPI */

#ifndef OLD_HDF5 /* prior to HDF5-1.6.5 */
#define hssize_t hsize_t
#endif

#define CLOUDY_METAL_COOLING 3

#endif
