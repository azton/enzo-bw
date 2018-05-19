/***********************************************************************
/
/  COMMUNICATION ROUTINE: RECEIVE FLUXES FROM ANOTHER PROCESSOR
/
/  written by: Greg Bryan
/  date:       December, 1997
/  modified1:
/
/  PURPOSE:
/
************************************************************************/
 
#include <stdio.h>
#include <string.h>
#include <assert.h>
 
#ifdef USE_MPI
#include "mpi.h"
#endif /* USE_MPI */
 
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "TopGridData.h"
#include "Hierarchy.h"
#include "LevelHierarchy.h"
 
 
 
 
int CommunicationReceiveFluxes(fluxes *Fluxes, int FromProc,
			       int NumberOfFields, int Rank)
{
 
  /* Count space and allocate buffer. */
 
  int dim1, dim2, field, i, TotalSize = 0, Sizes[MAX_DIMENSION], TempDim;
  for (dim1 = 0; dim1 < Rank; dim1++) {
    int size = 1;
    for (dim2 = 0; dim2 < Rank; dim2++) {
      TempDim = (Fluxes->LeftFluxEndGlobalIndex[dim1][dim2] -
	         Fluxes->LeftFluxStartGlobalIndex[dim1][dim2]) + 1;
      if (dim2 == dim1)
	TempDim = 1;
      size *= TempDim;
    }
    Sizes[dim1] = size;
    TotalSize += 2*size;
  }
 
  TotalSize *= NumberOfFields;
  float *buffer = new float[TotalSize];
 
  /* receive into buffer. */
 
#ifdef USE_MPI
 
  MPI_Status status;
  MPI_Datatype DataType = (sizeof(float) == 4) ? MPI_FLOAT : MPI_DOUBLE;
  MPI_Arg stat;
 
#ifdef MPI_INSTRUMENTATION
  starttime = MPI_Wtime();
#endif

  MPI_Arg Count = TotalSize;
  MPI_Arg Source = FromProc;

  if (MPI_Recv(buffer, Count, DataType, Source, MPI_FLUX_TAG, MPI_COMM_WORLD, &status) != MPI_SUCCESS) {
    fprintf(stderr, "Proc %"ISYM" MPI_Sendrecv error %"ISYM"\n", MyProcessorNumber, status.MPI_ERROR);
    return FAIL;
  }
 
#ifdef MPI_INSTRUMENTATION
  endtime = MPI_Wtime();
  timer[12] += endtime-starttime;
  counter[12] ++;
  RecvComm += endtime-starttime;
  CommunicationTime += endtime-starttime;
#endif /* MPI_INSTRUMENTATION */
 
 
#endif /* USE_MPI */
 
  /* Unpack buffer */
 
  int index = 0;
  for (dim1 = 0; dim1 < Rank; dim1++)
    for (field = 0; field < NumberOfFields; field++) {
      for (i = 0; i < Sizes[dim1]; i++)
	Fluxes->LeftFluxes[field][dim1][i] = buffer[index++];
      for (i = 0; i < Sizes[dim1]; i++)
	Fluxes->RightFluxes[field][dim1][i] = buffer[index++];
    }
 
  delete buffer;
 
  return SUCCESS;
}
