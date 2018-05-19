/***********************************************************************
/
/  GRID CLASS (DEPOSIT PARTICLE POSITIONS ONTO THE SPECIFIED FIELD)
/
/  written by: Greg Bryan
/  date:       May, 1995
/  modified1:
/
/  PURPOSE:
/     This routine deposits the particle living in this grid into either
/       the GravitatingMassField or the GravitatingMassFieldParticles of
/       the TargetGrid, depending on the value of DepositField.
/     It also moves the particles in this grid so they are at the correct
/       time and adjusts their 'mass' to be appropriate for TargetGrid
/       (since particle 'mass' changed from level to level).
/
/  NOTE:
/
************************************************************************/
 
#include <stdio.h>
#include <math.h>
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
 
/* This controls the maximum particle mass which will be deposited in
   the MASS_FLAGGING_FIELD.  Only set in Grid_SetFlaggingField. */
 
float DepositParticleMaximumParticleMass = 0;
 
int grid::DepositParticlePositions(grid *TargetGrid, FLOAT DepositTime,
				   int DepositField)
{
 
  /* Return if this doesn't concern us. */
 
  if (MyProcessorNumber != ProcessorNumber &&
      MyProcessorNumber != TargetGrid->ProcessorNumber)
    return SUCCESS;
 
  if (CommunicationDirection == COMMUNICATION_SEND &&
      (MyProcessorNumber != ProcessorNumber ||
       ProcessorNumber == TargetGrid->ProcessorNumber))
    return SUCCESS;
 
  if (CommunicationDirection == COMMUNICATION_RECEIVE &&
      MyProcessorNumber == ProcessorNumber &&
      ProcessorNumber != TargetGrid->ProcessorNumber)
    return SUCCESS;
 
//  fprintf(stderr, "----DPP: MyPN = %"ISYM", PN = %"ISYM", TGPN = %"ISYM", DIR (R=1,S=2) = %"ISYM", NP = %"ISYM"\n",
//    MyProcessorNumber, ProcessorNumber, TargetGrid->ProcessorNumber, CommunicationDirection, NumberOfParticles);
 
  /* Declarations. */
 
  int dim, i;
  float MassFactor = 1.0, *ParticleMassTemp, *ParticleMassPointer;
 
  /* If there are no particles, don't deposit anything. */
 
  if (NumberOfParticles == 0)
    return SUCCESS;
 
  /* If the particles are not on TargetGrid's processor, then transfer them. */
 
  if (ProcessorNumber != TargetGrid->ProcessorNumber) {
 
//  fprintf(stderr, "----DPP Call CommSendParticles PN = %"ISYM", TGPN = %"ISYM", NP = %"ISYM"\n",
//    ProcessorNumber, TargetGrid->ProcessorNumber, NumberOfParticles);
 
    if (this->CommunicationSendParticles(this, TargetGrid->ProcessorNumber, 0,
					 NumberOfParticles, 0) == FAIL) {
      fprintf(stderr, "Error in grid->CommunicationSendParticles.\n");
      return FAIL;
    }
    if (MyProcessorNumber == ProcessorNumber)
      return SUCCESS;   // Return, because we have finished transfering
  }
 
  /* If the Target is this grid and the DepositField is MassFlaggingField,
     then multiply the Particle density by the volume to get the mass. */
 
  if (this == TargetGrid && DepositField == MASS_FLAGGING_FIELD)
    for (dim = 0; dim < GridRank; dim++)
      MassFactor *= CellWidth[dim][0];
 
  /* If the DepositGrid and this grid are not the same, we must adjust the
     particle 'mass'. */
 
  if (this != TargetGrid) {
 
    /* Find the difference in resolution between this grid and TargetGrid. */
 
    float RefinementFactors[MAX_DIMENSION];
    this->ComputeRefinementFactorsFloat(TargetGrid, RefinementFactors);
 
    /* Compute the implied difference in 'mass' between particles in this
       grid and those in TargetGrid. */
 
    for (dim = 0; dim < GridRank; dim++)
      MassFactor *= RefinementFactors[dim];
 
  }
 
  /* If required, Change the mass of particles in this grid. */
 
  if (MassFactor != 1.0) {
    ParticleMassTemp = new float[NumberOfParticles];
    for (i = 0; i < NumberOfParticles; i++)
      ParticleMassTemp[i] = ParticleMass[i]*MassFactor;
    ParticleMassPointer = ParticleMassTemp;
  } else
    ParticleMassPointer = ParticleMass;
 
  /* If the target field is MASS_FLAGGING_FIELD, then set masses of
     particles which are too large to zero (to prevent run-away refinement). */
 
  if (DepositField == MASS_FLAGGING_FIELD &&
      DepositParticleMaximumParticleMass > 0 && MassFactor != 1.0)
    for (i = 0; i < NumberOfParticles; i++)
      ParticleMassPointer[i] = min(DepositParticleMaximumParticleMass,
				   ParticleMassPointer[i]);
 
  /* Compute difference between current time and DepositTime. */
 
  float TimeDifference = DepositTime - Time;
 
  /* Move particles to positions at Time + TimeDifference. */

  //The second argument forces the update even if
  //MyProcessor == Target->ProcessorNumber != this->ProcessorNumber
  this->UpdateParticlePosition(TimeDifference, TRUE);
 
  /* Deposit particles. */
 
//  fprintf(stderr, "----DPP Call TargetGrid->DepositPositions with NP = %"ISYM"\n", NumberOfParticles);
 
  if (TargetGrid->DepositPositions(ParticlePosition, ParticleMassPointer,
				   NumberOfParticles, DepositField) == FAIL) {
    fprintf(stderr, "Error in grid->DepositPositions\n");
    return FAIL;
  }
 
  /* If necessary, delete the particle mass temporary. */
 
  if (MassFactor != 1.0)
    delete ParticleMassTemp;
 
  /* If on a different processor, then we can just discard the particles. */
 
  if (ProcessorNumber != MyProcessorNumber) {
    this->DeleteAllFields();
    return SUCCESS;
  }
 
  /* Return particles to positions at Time. */
 
  this->UpdateParticlePosition(-TimeDifference);
 
  return SUCCESS;
}
