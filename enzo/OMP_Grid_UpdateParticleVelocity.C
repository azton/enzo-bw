/***********************************************************************
/
/  GRID CLASS (UPDATE PARTICLE VELOCITY FROM ACCELERATIONS)
/
/  written by: Greg Bryan
/  date:       May, 1995
/  modified1:  Robert Harkness
/  date:       9th August 2008
/              Simple OpenMP threads
/
/  PURPOSE:
/
/  NOTE:
/
************************************************************************/
 
#include <stdio.h>

#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
 
/* function prototypes */
 
#define VELOCITY_METHOD3
 
int CosmologyComputeExpansionFactor(FLOAT time, FLOAT *a, FLOAT *dadt);
 
 
int grid::UpdateParticleVelocity(float TimeStep)
{
 
  /* Return if this doesn't concern us. */
 
  if (ProcessorNumber != MyProcessorNumber)
    return SUCCESS;
 
  if (NumberOfParticles == 0 || SelfGravity == FALSE) return SUCCESS;
 
  FLOAT a = 1.0, dadt;
  int i;

#if defined(VELOCITY_METHOD1) || defined(VELOCITY_METHOD2)
  float VelocityMidStep;
#endif
 
  /* If using comoving coordinates, divide by a(t) first. */
 
  if (ComovingCoordinates)
    if (CosmologyComputeExpansionFactor(Time + TimeStep, &a, &dadt)
	== FAIL) {
      fprintf(stderr, "Error in CsomologyComputeExpansionFactors.\n");
      return FAIL;
    }
 
  /* Loop over dimensions. */
 
  for (int dim = 0; dim < GridRank; dim++) {
 
    /* Error check. */
 
    if (ParticleAcceleration[dim] == NULL) {
      fprintf(stderr, "No ParticleAccleration present.\n");
      return FAIL;
    }

  }
 
    /* Update velocities.  */
 
  if (ComovingCoordinates) {
 
    FLOAT coef = 0.5*dadt/a*TimeStep;
 
    /* If using comoving coordinates, subtract the (time-centered)
       drag-like term and add the acceleration. The acceleration has
        already been divided by a(t). */

    for (int dim = 0; dim < GridRank; dim++) {

#ifdef INTEL_OMP_SYNTAX
#pragma omp parallel private(i) shared(TimeStep, coef, dim) default(none)
#else
#pragma omp parallel private(i) shared(NumberOfParticles, ParticleVelocity, ParticleAcceleration, TimeStep, coef, dim) default(none)
#endif
    {

#pragma omp for schedule(static)
      for (i = 0; i < NumberOfParticles; i++) {
 
#ifdef VELOCITY_METHOD1
 
        /* i) partially time-centered. */
 
	VelocityMidStep = ParticleVelocity[dim][i] +
	                  ParticleAcceleration[dim][i]*0.5*TimeStep;
 
	ParticleVelocity[dim][i] +=
	  (-VelocityMidStep*dadt/a + ParticleAcceleration[dim][i]) * TimeStep;
 
#endif /* VELOCITY_METHOD1 */
 
#ifdef VELOCITY_METHOD2
 
        /* ii) partially backward. */
 
	VelocityMidStep = ParticleVelocity[dim][i] ;
 
	ParticleVelocity[dim][i] +=
	  (-VelocityMidStep*dadt/a + ParticleAcceleration[dim][i]) * TimeStep;
 
#endif /* VELOCITY_METHOD2 */
 
#ifdef VELOCITY_METHOD3
 
        /* iii) Semi-implicit way */
 
        ParticleVelocity[dim][i] = ((1.0 - coef)*ParticleVelocity[dim][i] +
                                    ParticleAcceleration[dim][i]*TimeStep)/
                                   (1.0 + coef);
 
#endif /* VELOCITY_METHOD3 */
 
      } // end loop over particles

    } // end parallel region

    } // end loop over dims

  } else {
 
      /* Otherwise, just add the acceleration. */

    for (int dim = 0; dim < GridRank; dim++) { 

#ifdef INTEL_OMP_SYNTAX
#pragma omp parallel private(i) shared(TimeStep, dim) default(none)
#else
#pragma omp parallel private(i) shared(NumberOfParticles, ParticleVelocity, ParticleAcceleration, TimeStep, dim) default(none)
#endif
    {

#pragma omp for schedule(static)
      for (i = 0; i < NumberOfParticles; i++) {
	ParticleVelocity[dim][i] += ParticleAcceleration[dim][i] * TimeStep;
      }

    } // end parallel region

    } // end loop over dims
 
  } // end if comoving etc.
 
  return SUCCESS;
}
