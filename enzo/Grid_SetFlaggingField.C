/***********************************************************************
/
/  GRID CLASS (CLEAR THE FLAGGING FIELD)
/
/  written by: Greg Bryan
/  date:       November, 1994
/  modified1:  Alexei Kritsuk, Aug. 2004: added refinement by shear.
/  modified2:  S.W. Skillman, Jan 2009:  Added refinement by shockwaves.
/
/  PURPOSE:
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
 
/* The following is defined in Grid_DepositParticlePositions.C. */
 
extern float DepositParticleMaximumParticleMass;
 
 
int grid::SetFlaggingField(int &NumberOfFlaggedCells, int level)
{
 
  /* Return if this doesn't concern us. */
 
  if (ProcessorNumber != MyProcessorNumber)
    return SUCCESS;
 
  /* declarations */
 
  NumberOfFlaggedCells = INT_UNDEFINED;
  int method;
 
  for (method = 0; method < MAX_FLAGGING_METHODS; method++) {
 
  /***********************************************************************/
  /* beginning of Cell flagging criterion routine                        */
 
    switch (CellFlaggingMethod[method]) {
 
    case 0:   /* no action */
      NumberOfFlaggedCells = (NumberOfFlaggedCells == INT_UNDEFINED ?
			      0 : NumberOfFlaggedCells);
      break;
 
      /* ==== METHOD 1: BY SLOPE ==== */
 
    case 1:
 
      /* flag all points needing extra resolution (FlagCellsToBeRefinedBySlop
	 returns the number of flagged cells). */
 
      NumberOfFlaggedCells = this->FlagCellsToBeRefinedBySlope();
      if (NumberOfFlaggedCells < 0) {
	fprintf(stderr, "Error in grid->FlagCellsToBeRefinedBySlope.\n");
	return FAIL;
      }
      break;
 
      /* ==== METHOD 2: BY BARYON MASS OR OVERDENSITY ==== */
 
    case 2:
 
      /* allocate and clear mass flagging field */
 
      this->ClearMassFlaggingField();
 
      /* baryons: add baryon density to mass flagging field (so the mass
	 flagging field contains the mass in the cell (not the density) */
 
      if (this->AddFieldMassToMassFlaggingField() == FAIL) {
	fprintf(stderr, "Error in grid->AddFieldMassToMassFlaggingField.\n");
	return FAIL;
      }
 
      /* flag all points that need extra resolution (FlagCellsToBeRefinedByMass
	 return the number of flagged cells). */
 
      NumberOfFlaggedCells = this->FlagCellsToBeRefinedByMass(level, method);
      if (NumberOfFlaggedCells < 0) {
	fprintf(stderr, "Error in grid->FlagCellsToBeRefinedByMass.\n");
	return FAIL;
      }
      break;
 
      /* ==== METHOD 3: BY SHOCKS ==== */
 
    case 3:
 
      NumberOfFlaggedCells = this->FlagCellsToBeRefinedByShocks();
      if (NumberOfFlaggedCells < 0) {
	fprintf(stderr, "Error in grid->FlagCellsToBeRefinedByShocks.\n");
	return FAIL;
      }
      break;
 
      /* ==== METHOD 4: BY PARTICLE MASS ==== */
 
    case 4:
 
      /* Allocate and clear mass flagging field. */
 
      this->ClearMassFlaggingField();
 
      /* Set the maximum particle mass to be deposited (cleared below). */
 
      DepositParticleMaximumParticleMass =
	0.99999*MinimumMassForRefinement[method]*POW(RefineBy,
		    level*MinimumMassForRefinementLevelExponent[method]);
 
      /* Deposit particles in this grid to MassFlaggingField. */
 
      if (this->DepositParticlePositions(this, this->ReturnTime(),
		  		         MASS_FLAGGING_FIELD) == FAIL) {
	fprintf(stderr, "Error in grid->DepositParticlePositions.\n");
	return FAIL;
      }
 
      DepositParticleMaximumParticleMass = 0;
 
      /* Flag all points that need extra resolution (FlagCellsToBeRefinedByMass
	 return the number of flagged cells). */
 
      NumberOfFlaggedCells = this->FlagCellsToBeRefinedByMass(level, method);
      if (NumberOfFlaggedCells < 0) {
	fprintf(stderr, "Error in grid->FlagCellsToBeRefinedByMass.\n");
	return FAIL;
      }
      break;
 
      /* ==== METHOD 6: BY JEANS LENGTH ==== */
 
    case 6:
 
      NumberOfFlaggedCells = this->FlagCellsToBeRefinedByJeansLength();
      if (NumberOfFlaggedCells < 0) {
	fprintf(stderr, "Error in grid->FlagCellsToBeRefinedByJeansLength.\n");
	return FAIL;
      }
      break;
 
      /* ==== METHOD 7: BY COOLING TIME < DX/SOUND SPEED ==== */
 
    case 7:
 
      NumberOfFlaggedCells = this->FlagCellsToBeRefinedByCoolingTime();
      if (NumberOfFlaggedCells < 0) {
	fprintf(stderr, "Error in grid->FlagCellsToBeRefinedByCoolingTime.\n");
	return FAIL;
      }
      break;
 
      /* ==== METHOD 8: BY POSITION OF MUST-REFINE PARTICLES  ==== */
 
    case 8:
 
      if (level < MustRefineParticlesRefineToLevel) {
	NumberOfFlaggedCells =
	  this->FlagCellsToBeRefinedByMustRefineParticles();
	if (NumberOfFlaggedCells < 0) {
	  fprintf(stderr, "Error in grid->FlagCellsByMustRefineParticles.\n");
	  return FAIL;
	}
      }
      if (NumberOfFlaggedCells == INT_UNDEFINED)
	NumberOfFlaggedCells = 0;
	
      break;
 
      /* ==== METHOD 9: BY SHEAR ==== */
 
    case 9:
 
      NumberOfFlaggedCells = this->FlagCellsToBeRefinedByShear();
      if (NumberOfFlaggedCells < 0) {
        fprintf(stderr, "Error in grid->FlagCellsToBeRefinedByShear.\n");
        return FAIL;
      }
      break;

      /* ==== METHOD 10: FORCE REFINEMENT TO SOME LEVEL IN A SET REGION ==== */
 
    case 10:
 
      NumberOfFlaggedCells = this->FlagCellsToBeRefinedByMustRefineRegion(level);
      if (NumberOfFlaggedCells < 0) {
        fprintf(stderr, "Error in grid->FlagCellsToBeRefinedByMustRefineRegion.\n");
        return FAIL;
      }
      break;
 

      /* ==== METHOD 11: FORCE REFINEMENT BASED ON METALLICITY OF GAS ==== */
 
    case 11:
 
      NumberOfFlaggedCells = this->FlagCellsToBeRefinedByMetallicity(level);
      if (NumberOfFlaggedCells < 0) {
        fprintf(stderr, "Error in grid->FlagCellsToBeRefinedByMetallicity.\n");
        return FAIL;
      }
      break;
 
//      /* ==== METHOD 12: Refine around Shockwaves (S.W. Skillman) ==== */
//
//   case 12:
//     NumberOfFlaggedCells = this->FlagCellsToBeRefinedByShockwaves();
//     if (NumberOfFlaggedCells < 0) {
//       fprintf(stderr, "Error in grid->FlagCellsToBeRefinedByShockwaves.\n");
//       return FAIL;
//     }
//     break;
      

#ifdef MHD
      /* ==== METHOD 20: BY HAND ==== */
      
    case 20:
      NumberOfFlaggedCells = this->FlagCellsToBeRefinedByHand(level);
      if (NumberOfFlaggedCells < 0) {
        fprintf(stderr, "Error in grid->FlagCellsToBeRefinedByHand.\n");
        return FAIL;
      }
#endif //MHD      


      /* ==== undefined ==== */
      
    case INT_UNDEFINED:
      break;
 
    default:
      fprintf(stderr, "CellFlaggingMethod[%"ISYM"] = %"ISYM" unknown\n", method,
	      CellFlaggingMethod[method]);
      return FAIL;
 
    }
 
    /* End of Cell flagging criterion routine                              */
    /***********************************************************************/
 
  } // end: loop over refinement methods
 
  if (NumberOfFlaggedCells == INT_UNDEFINED) {
    fprintf(stderr, "No valid CellFlaggingMethod specified.\n");
    return FAIL;
  }
 
#ifdef MPI_INSTRUMENTATION
  counter[4]++;
  timer[4] += NumberOfFlaggedCells;
#endif /* MPI_INSTRUMENTATION */
 
  if (debug1)
    printf("SetFlaggingField: NumberOfFlaggedCells = %"ISYM".\n",
	   NumberOfFlaggedCells);
 
  return SUCCESS;
 
}
