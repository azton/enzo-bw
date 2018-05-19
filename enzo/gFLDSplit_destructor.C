/*****************************************************************************
 *                                                                           *
 * Copyright 2009 Daniel R. Reynolds                                         *
 *                                                                           *
 * This software is released under the terms of the "Enzo Public License"    *
 * in the accompanying LICENSE file.                                         *
 *                                                                           *
 *****************************************************************************/
/***********************************************************************
/
/  Gray Flux-Limited Diffusion Split Implicit Problem Class 
/  Destructor routine
/
/  written by: Daniel Reynolds
/  date:       July 2009
/  modified1:  
/
/  PURPOSE: Frees all memory allocated for the implicit FLD problem.
/
************************************************************************/
#ifdef RAD_HYDRO
#include "gFLDSplit.h"



gFLDSplit::~gFLDSplit()
{

//   if (debug)  printf("Entering gFLDSplit::destructor routine\n");

  // delete HYPRE objects
  HYPRE_StructStencilDestroy(stencil);
  HYPRE_StructGridDestroy(grid);

  // delete EnzoVectors and other internal arrays
  int i, j;
  //   EnzoVectors require deleting the structure
  delete sol;
  delete U0;
  delete extsrc;

  //   arrays require deleting the array
  HYPRE_StructVectorDestroy(rhsvec);
  HYPRE_StructVectorDestroy(solvec);
  HYPRE_StructMatrixDestroy(P);
  delete[] matentries;
  delete[] rhsentries;
  delete[] HYPREbuff;
  delete[] FluidEnergyCorrection;
  delete[] OpacityE;
  delete[] Temperature;
  delete[] Temperature0;

  // delete boundary condition arrays
  for (i=0; i<3; i++)
    for (j=0; j<2; j++) {
      if (BdryVals[i][j] != NULL)  delete[] BdryVals[i][j];
    }

}
#endif
