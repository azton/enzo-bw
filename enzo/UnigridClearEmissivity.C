/*********************************************************************
/
/  CLEAR EMISSIVITY FUNCTION
/
/  written by: Geoffrey So
/  date:       Jan 9, 2009
/              Jan 15, 2010 cleaned it up to match ClearEmissivity.C
/        
/  PURPOSE:
/    To clear the emissivity array for Unigrid = 0
/    (implimented as a BaryonField) when
/    starting in the next timestep in the root grid, so that when
/    CalcEmiss(...) is called, BaryonField[EmisNum] will have all values 
/    initalized to 0
/
/
/
/
*********************************************************************/

#include <stdio.h>
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "fortran.def"

int FindField(int field, int farray[], int numfields);

#ifdef EMISSIVITY
int grid::UnigridClearEmissivity(){
  if( MyProcessorNumber != ProcessorNumber)
    return SUCCESS;
  int size = GridDimension[0]*GridDimension[1]*GridDimension[2];
  int EtaNum;

  // Copying what's done in Grid_IdentifyPhysicalQuantities to check
  if ((EtaNum = FindField(EmissivityField0, FieldType, NumberOfBaryonFields)) < 0) {
    fprintf(stderr, "Cannot find EmissivityField0.\n");
    return FAIL;
  }

  //EtaNum = FindField(EmissivityField0, FieldType, NumberOfBaryonFields);
  for(int i=0; i<size; i++)
    BaryonField[EtaNum][i] = 0;
  //  printf("DONE \n");
  return SUCCESS;
}
#endif
