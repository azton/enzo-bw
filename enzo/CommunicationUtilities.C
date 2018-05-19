/***********************************************************************
/
/  COMMUNICATION ROUTINE: FIND MINIMUM VALUE AMOUNG PROCESSORS
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
 
 
 
 
float CommunicationMinValue(float Value)
{
 
  if (NumberOfProcessors == 1)
    return Value;
 
  float ReturnValue = Value;
 
#ifdef USE_MPI
 
  MPI_Datatype DataType = (sizeof(float) == 4) ? MPI_FLOAT : MPI_DOUBLE;
  MPI_Arg Count;
 
//  printf("min: %"ISYM" sending %"FSYM"\n", MyProcessorNumber, Value);
 
#ifdef MPI_INSTRUMENTATION
  starttime = MPI_Wtime();
#endif

  Count = 1;
 
  MPI_Allreduce(&Value, &ReturnValue, Count, DataType, MPI_MIN, MPI_COMM_WORLD);
 
#ifdef MPI_INSTRUMENTATION
  endtime = MPI_Wtime();
  timer[16]+= endtime-starttime;
  counter[16] ++;
  GlobalCommunication += endtime-starttime;
  CommunicationTime += endtime-starttime;
#endif /* MPI_INSTRUMENTATION */
 
#endif /* USE_MPI */
 
  return ReturnValue;
}
 
 
float CommunicationMaxValue(float Value)
{
 
  if (NumberOfProcessors == 1)
    return Value;
 
  float ReturnValue = Value;
 
#ifdef USE_MPI
 
  MPI_Datatype DataType = (sizeof(float) == 4) ? MPI_FLOAT : MPI_DOUBLE;
  MPI_Arg Count;
 
//  printf("min: %"ISYM" sending %"FSYM"\n", MyProcessorNumber, Value);
 
#ifdef MPI_INSTRUMENTATION
  starttime = MPI_Wtime();
#endif

  Count = 1;
 
  MPI_Allreduce(&Value, &ReturnValue, Count, DataType, MPI_MAX, MPI_COMM_WORLD);
 
#ifdef MPI_INSTRUMENTATION
  endtime = MPI_Wtime();
  timer[16]+= endtime-starttime;
  counter[16] ++;
  GlobalCommunication += endtime-starttime;
  CommunicationTime += endtime-starttime;
#endif /* MPI_INSTRUMENTATION */
 
#endif /* USE_MPI */
 
  return ReturnValue;
}
 
 
int CommunicationSumValues(float *Values, int Number)
{
 
  if (NumberOfProcessors == 1)
    return SUCCESS;
 
#ifdef USE_MPI
 
  MPI_Datatype DataType = (sizeof(float) == 4) ? MPI_FLOAT : MPI_DOUBLE;
  MPI_Arg Count;
 
  int i;
  float *buffer = new float[Number];
  for (i = 0; i < Number; i++)
    buffer[i] = Values[i];
 
  double starttime = MPI_Wtime();

  Count = Number;

  MPI_Reduce(buffer, Values, Count, DataType, MPI_SUM, ROOT_PROCESSOR, MPI_COMM_WORLD);
 
  double endtime = MPI_Wtime();
 
  delete [] buffer;
 
  CommunicationTime += endtime-starttime;
 
#endif /* USE_MPI */
 
  return SUCCESS;
}
 
 
 
 
int CommunicationAllSumValues(float *Values, int Number)
{
 
  if (NumberOfProcessors == 1)
    return SUCCESS;
 
#ifdef USE_MPI
 
  MPI_Datatype DataType = (sizeof(float) == 4) ? MPI_FLOAT : MPI_DOUBLE;
  MPI_Arg Count;
 
  int i;
  float *buffer = new float[Number];
  for (i = 0; i < Number; i++)
    buffer[i] = Values[i];
 
  double starttime = MPI_Wtime();

  Count = Number;
 
  MPI_Allreduce(buffer, Values, Count, DataType, MPI_SUM, MPI_COMM_WORLD);
 
  double endtime = MPI_Wtime();
 
  delete [] buffer;
 
  CommunicationTime += endtime-starttime;
 
#endif /* USE_MPI */
 
  return SUCCESS;
}
 
 
 
 
int CommunicationAllSumIntegerValues(int *Values, int Number)
{
  if (NumberOfProcessors == 1)
    return SUCCESS;
 
#ifdef USE_MPI

  MPI_Datatype DataTypeInt = (sizeof(int) == 4) ? MPI_INT : MPI_LONG_LONG_INT; 
  MPI_Arg Count;
 
  int i;
  int *buffer = new int[Number];
  for (i = 0; i < Number; i++)
    buffer[i] = Values[i];
 
  double starttime = MPI_Wtime();

  Count = Number;
 
  MPI_Allreduce(buffer, Values, Count, DataTypeInt, MPI_SUM, MPI_COMM_WORLD);
 
  double endtime = MPI_Wtime();
 
  delete [] buffer;
 
  CommunicationTime += endtime-starttime;
 
#endif /* USE_MPI */
 
  return SUCCESS;
}
