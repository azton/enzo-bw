/***********************************************************************
/
/  GRID CLASS (DEPOSIT BARYON FIELD IN TO TARGET GRAVITATING MASS FIELD)
/
/  written by: Greg Bryan
/  date:       March, 1999
/  modified1:  Robert Harkness
/  date:       March, 2004
/
/  PURPOSE:
/
/  RETURNS: FAIL or SUCCESS
/
************************************************************************/
 
#include <stdio.h>
 
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
 
int CosmologyComputeExpansionFactor(FLOAT time, FLOAT *a, FLOAT *dadt);
 
#ifdef USE_MPI
int CommunicationBufferedSend(void *buffer, int size, MPI_Datatype Type, int Target,
			      int Tag, MPI_Comm CommWorld, int BufferSize);
#endif /* USE_MPI */
 
extern "C" void FORTRAN_NAME(dep_grid_cic)(
                               float *source, float *dest, float *temp,
			       float *velx, float *vely, float *velz,
			       float *dt, float *rfield, int *ndim,
                                   hydro_method *ihydro,
			       float *delx, float *dely, float *delz,
			       int *sdim1, int *sdim2, int *sdim3,
			       int *sstart1, int *sstart2, int *sstart3,
			       int *send1, int *send2, int *send3,
			       float *offset1, float *offset2, float *offset3,
			       int *ddim1, int *ddim2, int *ddim3,
			       int *refine1, int *refine2, int *refine3);
 
/* InterpolateBoundaryFromParent function */
 
int grid::DepositBaryons(grid *TargetGrid, FLOAT DepositTime)
{
 
  /* If this doesn't concern us, return. */
 
  if ((MyProcessorNumber != ProcessorNumber &&
       MyProcessorNumber != TargetGrid->ProcessorNumber) ||
       NumberOfBaryonFields == 0)
    return SUCCESS;
 
  if (CommunicationDirection == COMMUNICATION_SEND &&
      (MyProcessorNumber != ProcessorNumber ||
       ProcessorNumber == TargetGrid->ProcessorNumber))
    return SUCCESS;
 
  if (CommunicationDirection == COMMUNICATION_RECEIVE &&
      MyProcessorNumber == ProcessorNumber &&
      ProcessorNumber != TargetGrid->ProcessorNumber)
    return SUCCESS;
 
  TargetGrid->DebugCheck("DepositBaryons_target");
  this->DebugCheck("DepositBaryons_this");
 
  /* Declarations. */
 
  float GridStart[MAX_DIMENSION] = {0,0,0};
  int GridOffset[MAX_DIMENSION] = {0,0,0}, Refinement[MAX_DIMENSION],
      RegionDim[MAX_DIMENSION] = {1,1,1}, GridOffsetEnd[MAX_DIMENSION],
      i, j, k, index, gmindex, dim, size = 1;
 
  /* Find fields: density, total energy, velocity1-3. */
 
  int DensNum, GENum, Vel1Num, Vel2Num, Vel3Num, TENum;
  if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num,
					 Vel3Num, TENum) == FAIL) {
    fprintf(stderr, "Error in IdentifyPhysicalQuantities.\n");
    return FAIL;
  }
 
  /* Error check: subgrid covering field must exist on entry. */
 
  if (MyProcessorNumber == ProcessorNumber &&
      BaryonField[NumberOfBaryonFields] == NULL) {
    fprintf(stderr, "subgrid covering field missing\n");
    return FAIL;
  }
 
  /* Compute refinement factors. */
 
  TargetGrid->ComputeRefinementFactors(this, Refinement);
 
  /* This routine will create a temporary patch with cell width equal to
     the target grid.  The current grid then deposits into this patch.
     Compute the TargetOffset (in grid units) and TargetStartIndex and
     the region dim (in Target units). */
 
  for (dim = 0; dim < GridRank; dim++) {
 
    /* GridOffset is the number of TargetGrid cells from the edge of
       the TargetGrid mass field and the region to be deposited by this
       grid.  It must not extended beyond the active region of TargetGrid
       (if we are depositing in self). */
 
    GridOffset[dim] = nint((GridLeftEdge[dim] -
			    TargetGrid->GravitatingMassFieldLeftEdge[dim])/
			   TargetGrid->GravitatingMassFieldCellSize) - 1;
    if (TargetGrid == this)
      GridOffset[dim] = max(GridOffset[dim],
	nint((TargetGrid->GridLeftEdge[dim] -
	      TargetGrid->GravitatingMassFieldLeftEdge[dim])/
	     TargetGrid->GravitatingMassFieldCellSize)   );
 
    /* GridStart is the distance (float) in target grid cell units between
       the exact value of GridOffset and it's integer version. */
 
    GridStart[dim] = (GridLeftEdge[dim] -
		      TargetGrid->GravitatingMassFieldLeftEdge[dim])/
                      TargetGrid->GravitatingMassFieldCellSize -
                      float(GridOffset[dim]);
 
    /* RegionDim is the size, in TargetGrid cell units, of the region (patch)
       to be deposited. It must not extend beyond the edge of the active
       region of TargetGrid. */
 
    GridOffsetEnd[dim] = nint((GridRightEdge[dim] -
			    TargetGrid->GravitatingMassFieldLeftEdge[dim])/
			   TargetGrid->GravitatingMassFieldCellSize);
    if (TargetGrid == this)
      GridOffsetEnd[dim] = min(GridOffsetEnd[dim],
	nint((TargetGrid->GridRightEdge[dim] -
	      TargetGrid->GravitatingMassFieldLeftEdge[dim])/
	     TargetGrid->GravitatingMassFieldCellSize)-1 );
    RegionDim[dim] = GridOffsetEnd[dim] - GridOffset[dim] + 1;
		
    size *= RegionDim[dim];
 
    if (TargetGrid != this && GridStart[dim] < 0) {
      fprintf(stderr, "GridStart[%"ISYM"] = %"GSYM" < 0.\n", dim,GridStart[dim]);
      return FAIL;
    }
 
    if (RegionDim[dim] < 2) {
      fprintf(stderr, "RegionDim[%"ISYM"] = %"ISYM" < 2\n", dim, RegionDim[dim]);
      return FAIL;
    }
 
  }
 
  /* Prepare the density field. */
 
  float *dens_field = new float[size];
 
  if (ProcessorNumber == MyProcessorNumber) {
 
    /* Compute the dt to advance from current time to DepositTime. */
 
    FLOAT a = 1, dadt;
    if (ComovingCoordinates)
      if (CosmologyComputeExpansionFactor(0.5*(Time+DepositTime), &a, &dadt)
	  == FAIL) {
      fprintf(stderr, "Error in CosmologyComputeExpansionFactors.\n");
      return FAIL;
    }
    float dt = (DepositTime - Time)/a;
    dt = 0;
 
    /* Set up a float version of cell size to pass to fortran. */
 
    float dxfloat[MAX_DIMENSION] = {0,0,0};
    for (dim = 0; dim < GridRank; dim++)
      dxfloat[dim] = float(CellWidth[dim][0]);
 
    /* Allocate a density and velocity mesh for this grid. */
 
    float *vel_field = new float[size*4];
 
    /* Generate the density field advanced by dt using smoothed
       velocity field. */
 
//  fprintf(stderr, "Grid_DepositBaryons - call dep_grid_cic\n");
 
    FORTRAN_NAME(dep_grid_cic)(BaryonField[DensNum], dens_field, vel_field,
			 BaryonField[Vel1Num], BaryonField[Vel2Num],
			     BaryonField[Vel3Num], &dt,
			 BaryonField[NumberOfBaryonFields], &GridRank,
                             &HydroMethod,
			 dxfloat, dxfloat+1, dxfloat+2,
			 GridDimension, GridDimension+1, GridDimension+2,
			 GridStartIndex, GridStartIndex+1, GridStartIndex+2,
			 GridEndIndex, GridEndIndex+1, GridEndIndex+2,
			 GridStart, GridStart+1, GridStart+2,
			 RegionDim, RegionDim+1, RegionDim+2,
			 Refinement, Refinement+1, Refinement+2);
 
    delete [] vel_field;
 
  } // end: if (ProcessorNumber == MyProcessorNumber)
 
  /* If necessary, copy data from this processor to target grid's processor.
     Note: this really needs to be put into it's own transfer routine. */
 
  if (ProcessorNumber != TargetGrid->ProcessorNumber) {
 
#ifdef USE_MPI
 
    MPI_Request  RequestHandle;
    MPI_Status status;
    MPI_Datatype DataType = (sizeof(float) == 4) ? MPI_FLOAT : MPI_DOUBLE;
    MPI_Arg Count;
    MPI_Arg Source;

    Count = size;
    Source = ProcessorNumber;
 
    double time1 = MPI_Wtime();
 
    if (MyProcessorNumber == TargetGrid->ProcessorNumber)
      MPI_Irecv(dens_field, Count, DataType, Source,
               MPI_SENDREGION_TAG, MPI_COMM_WORLD, &RequestHandle);
 
 
    if (MyProcessorNumber == ProcessorNumber)
      CommunicationBufferedSend(dens_field, size, DataType,
			     TargetGrid->ProcessorNumber, MPI_SENDREGION_TAG,
				MPI_COMM_WORLD, BUFFER_IN_PLACE);
 
    double time2 = MPI_Wtime();
 
    if (MyProcessorNumber == TargetGrid->ProcessorNumber)
      MPI_Wait(&RequestHandle, &status);
 
    double time3 = MPI_Wtime();
 
    CommunicationTime += time3 - time1;
    WaitComm += time3 - time2;
 
#endif /* USE_MPI */
 
  } // end: if (ProcessorNumber != TargetGrid->ProcessorNumber)
 
  /* Return if this is not our concern. */
 
  if (MyProcessorNumber != TargetGrid->ProcessorNumber) {
    //    delete [] dens_field;  new done in CommBufferedSend
    delete [] BaryonField[NumberOfBaryonFields];
    BaryonField[NumberOfBaryonFields] = NULL;
    return SUCCESS;
  }
 
  /* Add dens_field to GravitatingMassField in target grid. */
 
  index = 0;
  for (k = 0; k < RegionDim[2]; k++)
    for (j = 0; j < RegionDim[1]; j++) {
      gmindex = (j+GridOffset[1] +
               (k+GridOffset[2])*TargetGrid->GravitatingMassFieldDimension[1])*
	      TargetGrid->GravitatingMassFieldDimension[0] + GridOffset[0];
      for (i = 0; i < RegionDim[0]; i++, gmindex++, index++)
	TargetGrid->GravitatingMassField[gmindex] += dens_field[index];
    }
 
  /* Clean up (first under_subgrid field). */
 
  delete [] dens_field;
  delete [] BaryonField[NumberOfBaryonFields];
  BaryonField[NumberOfBaryonFields] = NULL;
 
  return SUCCESS;
 
}
