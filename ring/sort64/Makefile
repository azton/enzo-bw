#=======================================================================
#  MAKEFILE FOR RING I/O
#  Robert Harkness, 26th June 2006
#=======================================================================
#
# ENZO_INT is SMALL_INTS for 32 bits or LARGE_INTS for 64 bits
# ENZO_XPR is r4 (32 bits) or r8 (64 bits) for float and real

null :=
hdf5 := hdf5

ENZO_PPR = p8
ENZO_XPR = r8
ENZO_INT = LARGE_INTS
#ENZO_INT = SMALL_INTS

SYSTEM := $(shell uname -s)
UNAMEM := $(shell uname -m)

ifeq ($(ENZO_ARCH),$(null))

ifeq ($(SYSTEM),AIX)
  ifeq ($(OBJECT_MODE),$(null))
    ENZO_ARCH = aix32
    ENZO_MESS = mpi
  endif
  ifeq ($(OBJECT_MODE),64)
    ENZO_ARCH = aix64
    ENZO_MESS = mpi
  endif
endif

ifeq ($(SYSTEM),Linux)
  ifeq ($(UNAMEM),x86_64)
    ENZO_ARCH = xt5
    ENZO_MESS = mpi
  endif
endif



ifeq ($(SYSTEM),Darwin) 
  ifeq ($(UNAMEM),x86_64)
    ENZO_ARCH = geoffrey
    ENZO_MESS = mpi
  endif
  ifeq ($(ENZO_INT),LARGE_INTS)
    ifeq ($(ENZO_XPR),r8)
      FORTPR = -fdefault-real-8 -fdefault-integer-8
    else
      FORTPR = -fdefault-integer-8
    endif
  else
    ifeq ($(ENZO_XPR),r8)
      FORTPR = -fdefault-real-8
    else
      FORTPR =
    endif
  endif  
endif









endif


ifeq ($(ENZO_HDF),$(null))
  ENZO_HDF := hdf5
endif

ifeq ($(ENZO_END),$(null))
  ENZO_END := HDF5_BE
endif

ifeq ($(SYSTEM),Linux)
  CPP_FLAG := -traditional
endif

#-------------------  object files  ------------------------------------

OBJ = PF_Ring_Decomp.o \
	Enzo_Dims_create.o \
	Mpich_V1_Dims_create.o

ifeq ($(ENZO_INT),SMALL_INTS)
  EXENAME := Sort_32_PF
else
  EXENAME := Sort_64_PF
endif


#-------------------  macro definitions -------------------------------

.SUFFIXES:
.SUFFIXES: .src .c .f .o .C

HEADER  = macros_and_parameters.h extern_hdf5.h
EXEDIR  =

#------------------ include machine-dependent macros -----------------

INCLUDE_MAKEFILE = Make.${ENZO_ARCH}.${ENZO_MESS}
include ${INCLUDE_MAKEFILE}

#--------------------  implicit rules ---------------------------------

.src.o:
	${CPP} ${DEFINES} -P ${CPP_FLAG} $< > $*.f
	${FC} ${FC_OPT} ${DEBUG} ${OPTS} $*.f 
.C.o:
	${CCxx} ${CCxx_OPT} -D${ENZO_INT} -D${ENZO_XPR} -D${ENZO_END} ${DEBUG} ${OPTS} $*.C
.c.o:
	${CC} ${CC_OPT} -D${ENZO_INT} -D${ENZO_XPR} -D${ENZO_END} ${DEBUG} ${OPTS} -c $<

#-------------------------  main targets ------------------------------

${EXENAME} : ${OBJ}
	${LDR} ${LDR_OPT} ${DEBUG} -o ${EXEDIR}${EXENAME} ${OBJ} ${LIB}

#---------------------  special compilation rules ---------------------


#----------------------------- help ------------------------------------

help:
	@echo Type 'make' to generate executable
	@echo Type 'make clean' to remove object files, etc

#----------------------------- clean ------------------------------------

clean:
	rm -f *.o *.f ${EXEDIR}${EXENAME}

#-----------------------------------------------------------------------

${OBJ}: ${HEADER}
