/***********************************************************************
/
/  GRID CLASS (SEND PARTICLES FROM REAL GRID TO 'FAKE' (REPLICATED) GRID)
/
/  written by: Greg Bryan
/  date:       December, 1997
/  modified1:  Robert Harkness
/  date:       March, 2004
/
/  PURPOSE:
/
/  INPUTS:
/
************************************************************************/
 
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
 
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
 
/* function prototypes */
 
#ifdef USE_MPI
int CommunicationBufferedSend(void *buffer, int size, MPI_Datatype Type, int Target,
			      int Tag, MPI_Comm CommWorld, int BufferSize);
#endif /* USE_MPI */
 
/* Send particle from this grid to ToGrid on processor ToProcessor, using
   FromNumber particles counting from FromStart.  Place into ToGrid at
   particle number ToStart. If ToStart = -1, then add to end. */
 
int grid::CommunicationSendParticles(grid *ToGrid, int ToProcessor,
				    int FromStart, int FromNumber, int ToStart)
{
 
  int i, j, dim, index;
 
  if (FromNumber == 0)
    return SUCCESS;
 
  /* Compute size of region to transfer. */
 
  int NumberOfFields = 1+GridRank+NumberOfParticleAttributes;
  int RegionSize = FromNumber;
  int TransferSize = RegionSize * NumberOfFields;
 
  //  fprintf(stderr, "P(%"ISYM"): sending %"ISYM" from %"ISYM" -> %"ISYM" (%"ISYM" %"ISYM")\n",
  //  	  MyProcessorNumber, FromNumber, ProcessorNumber, ToProcessor,
  //  	  FromStart, ToStart);
 
  /* Allocate buffer. */
 
  float *buffer = NULL;
  if (MyProcessorNumber == ProcessorNumber ||
      MyProcessorNumber == ToProcessor)
    buffer = new float[TransferSize];
 
  /* If this is the from processor, pack fields. */
 
  if (MyProcessorNumber == ProcessorNumber) {
 
    index = 0;
    for (i = FromStart; i < FromStart+FromNumber; i++)
      buffer[index++] = ParticleMass[i];
    for (dim = 0; dim < GridRank; dim++)
      for (i = FromStart; i < FromStart+FromNumber; i++) {
	//	buffer[index++] = ParticlePosition[dim][i];
	buffer[index++] = ParticleVelocity[dim][i];
      }
    for (j = 0; j < NumberOfParticleAttributes; j++)
      for (i = FromStart; i < FromStart+FromNumber; i++)
	buffer[index++] = ParticleAttribute[j][i];
 
    if (index != TransferSize) {
      fprintf(stderr, "index = %"ISYM"  TransferSize = %"ISYM"\n", index, TransferSize);
      return FAIL;
    }
 
  } // end: if (MyProcessorNumber)
 
  /* Allocate Number field on from processor. */
 
  FLOAT *TempPos[MAX_DIMENSION];
  float  *TempVel[MAX_DIMENSION], *TempMass,
        *TempAttribute[MAX_NUMBER_OF_PARTICLE_ATTRIBUTES];
  int *TempNumber, NewNumber = FromNumber, *TempType;
  if (ToStart == -1)
    NewNumber += ToGrid->NumberOfParticles;
 
  if (MyProcessorNumber == ToProcessor) {
 
    /* If ToStart == -1, then add new particles to end (first, copy to buffer,
       then allocate and copy back. */
 
    if (ToStart == -1) {
      TempMass = ToGrid->ParticleMass;
      TempNumber = ToGrid->ParticleNumber;
      TempType = ToGrid->ParticleType;
      for (dim = 0; dim < MAX_DIMENSION; dim++) {
	TempPos[dim] = ToGrid->ParticlePosition[dim];
	TempVel[dim] = ToGrid->ParticleVelocity[dim];
      }
      for (j = 0; j < NumberOfParticleAttributes; j++)
	TempAttribute[j] = ToGrid->ParticleAttribute[j];
      ToGrid->ParticleNumber = NULL;  // signal that we should reallocate
    }
 
    /* If unallocated, then allocate. */
 
    if (ToGrid->ParticleNumber == NULL) {
      ToGrid->ParticleNumber = new int[NewNumber];
      ToGrid->ParticleMass   = new float[NewNumber];
      ToGrid->ParticleType   = new int[NewNumber];
      for (dim = 0; dim < MAX_DIMENSION; dim++) {
	ToGrid->ParticlePosition[dim] = new FLOAT[NewNumber];
	ToGrid->ParticleVelocity[dim] = new float[NewNumber];
      }
      for (j = 0; j < NumberOfParticleAttributes; j++)
	ToGrid->ParticleAttribute[j] = new float[NewNumber];
      if (ToStart > 0) {
	fprintf(stderr, "Unallocated Number, yet FromStart = %"ISYM"\n", FromStart);
	return FAIL;
      }
    }
 
    /* If adding to end, then copy and delete old fields. */
 
    if (ToStart == -1) {
      for (i = 0; i < ToGrid->NumberOfParticles; i++) {
	ToGrid->ParticleNumber[i] = TempNumber[i];
	ToGrid->ParticleMass[i]   = TempMass[i];
	ToGrid->ParticleType[i]   = TempType[i];
      }
      for (dim = 0; dim < GridRank; dim++)
	for (i = 0; i < ToGrid->NumberOfParticles; i++) {
	  ToGrid->ParticlePosition[dim][i] = TempPos[dim][i];
	  ToGrid->ParticleVelocity[dim][i] = TempVel[dim][i];
	}
      for (j = 0; j < NumberOfParticleAttributes; j++)
	for (i = 0; i < ToGrid->NumberOfParticles; i++)
	  ToGrid->ParticleAttribute[j][i] = TempAttribute[j][i];
	
      delete [] TempNumber;
      delete [] TempMass;
      delete [] TempType;
      for (dim = 0; dim < GridRank; dim++) {
	delete [] TempPos[dim];
	delete [] TempVel[dim];
      }
      for (j = 0; j < NumberOfParticleAttributes; j++)
	delete [] TempAttribute[j];
      ToStart = ToGrid->NumberOfParticles;
    }
 
  } // end: if (MyProcessorNumber == ToProcessor)
 
  ToGrid->NumberOfParticles = max(NewNumber, ToGrid->NumberOfParticles);
 
  /* Send buffer. */
 
#ifdef USE_MPI
 
  /* only send if processor numbers are not identical */
 
  if (ProcessorNumber != ToProcessor) {
 
    MPI_Request RequestHandle[6];
    MPI_Status status;
    MPI_Datatype DataType = (sizeof(float) == 4) ? MPI_FLOAT : MPI_DOUBLE;
    MPI_Datatype DataTypeInt = (sizeof(int) == 4) ? MPI_INT : MPI_LONG_LONG_INT;
    MPI_Arg Count;
    MPI_Arg Source;

 
    int irq;
 
#ifdef MPI_INSTRUMENTATION
    starttime = MPI_Wtime();
#endif
 
    if (MyProcessorNumber == ToProcessor) {
      irq = 0;

      Count = TransferSize;
      Source = ProcessorNumber;
 
      MPI_Irecv(buffer, Count, DataType, Source,
                MPI_SENDPARTFLOAT_TAG, MPI_COMM_WORLD, &RequestHandle[irq++]);

      Count = FromNumber;
      Source = ProcessorNumber;
 
      MPI_Irecv(&ToGrid->ParticleNumber[ToStart], Count, DataTypeInt, Source,
                MPI_SENDPARTINT_TAG, MPI_COMM_WORLD, &RequestHandle[irq++]);

      Count = FromNumber;
      Source = ProcessorNumber;
 
      MPI_Irecv(&ToGrid->ParticleType[ToStart], Count, DataTypeInt, Source,
                MPI_SENDPARTINT_TAG, MPI_COMM_WORLD, &RequestHandle[irq++]);
 
      for (dim = 0; dim < GridRank; dim++) {
        Count = FromNumber;
        Source = ProcessorNumber;

        MPI_Irecv(ToGrid->ParticlePosition[dim]+ToStart, Count, MY_MPIFLOAT, Source,
                  MPI_SENDPARTFLOAT_TAG, MPI_COMM_WORLD, &RequestHandle[irq++]);
      }

    }
 
    if (MyProcessorNumber == ProcessorNumber) {
 
      CommunicationBufferedSend(buffer, TransferSize, DataType, ToProcessor,
	       MPI_SENDPARTFLOAT_TAG, MPI_COMM_WORLD, BUFFER_IN_PLACE);
 
      CommunicationBufferedSend(ParticleNumber+FromStart, FromNumber, DataTypeInt,
				ToProcessor, MPI_SENDPARTINT_TAG,
				MPI_COMM_WORLD, FromNumber*sizeof(int));
 
      CommunicationBufferedSend(ParticleType+FromStart, FromNumber, DataTypeInt,
				ToProcessor, MPI_SENDPARTINT_TAG,
				MPI_COMM_WORLD, FromNumber*sizeof(int));
 
      for (dim = 0; dim < GridRank; dim++)
	CommunicationBufferedSend(ParticlePosition[dim]+FromStart, FromNumber, MY_MPIFLOAT,
				  ToProcessor, MPI_SENDPARTFLOAT_TAG,
			          MPI_COMM_WORLD, FromNumber*sizeof(FLOAT));
    }
 
#ifdef MPI_INSTRUMENTATION
    endtime = MPI_Wtime();
    timer[7] += endtime-starttime;
    counter[7] ++;
    RecvComm += endtime - starttime;
    CommunicationTime += endtime - starttime;
#endif /* MPI_INSTRUMENTATION */
 
#ifdef MPI_INSTRUMENTATION
    starttime = MPI_Wtime();
#endif
 
    if (MyProcessorNumber == ToProcessor) {
      irq = 0;
 
      MPI_Wait(&RequestHandle[irq++], &status);
 
      MPI_Wait(&RequestHandle[irq++], &status);
 
      MPI_Wait(&RequestHandle[irq++], &status);
 
      for (dim = 0; dim < GridRank; dim++)
	MPI_Wait(&RequestHandle[irq++], &status);
    }
 
#ifdef MPI_INSTRUMENTATION
    endtime = MPI_Wtime();
    timer[7] += endtime-starttime;
    WaitComm += endtime - starttime;
    CommunicationTime += endtime - starttime;
#endif /* MPI_INSTRUMENTATION */
 
 
  } // end: if (ProcessorNumber != ToProcessor)
 
#endif /* USE_MPI */
 
  /* If this is the to processor, unpack fields. */
 
  if (MyProcessorNumber == ToProcessor) {
 
    index = 0;
    for (i = ToStart; i < ToStart+FromNumber; i++)
      ToGrid->ParticleMass[i] = buffer[index++];
    for (dim = 0; dim < GridRank; dim++)
      for (i = ToStart; i < ToStart+FromNumber; i++) {
	//	ToGrid->ParticlePosition[dim][i] = buffer[index++];
	ToGrid->ParticleVelocity[dim][i] = buffer[index++];
      }
    for (j = 0; j < NumberOfParticleAttributes; j++)
      for (i = ToStart; i < ToStart+FromNumber; i++)
	ToGrid->ParticleAttribute[j][i] = buffer[index++];
 
    /* If on the same processor, don't forget to copy the number field. */
 
    if (ProcessorNumber == ToProcessor && this != ToGrid) {
      for (i = 0; i < FromNumber; i++) {
	ToGrid->ParticleNumber[i+ToStart] = ParticleNumber[i+FromStart];
        ToGrid->ParticleType[i+ToStart] = ParticleType[i+FromStart];
      }
      for (dim = 0; dim < GridRank; dim++)
        for (i = 0; i < FromNumber; i++)
          ToGrid->ParticlePosition[dim][i+ToStart] = ParticlePosition[dim][i+FromStart];
    }
			
  } // end: if (MyProcessorNumber...)
 
  if (MyProcessorNumber == ToProcessor)
    delete [] buffer;
 
  return SUCCESS;
}
 
