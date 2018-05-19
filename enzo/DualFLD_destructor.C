/*****************************************************************************
 * Copyright 2013 Daniel R. Reynolds                                         *
 *                                                                           *
 * This software is released under the terms of the "Enzo Public License"    *
 * in the accompanying LICENSE file.                                         *
 *****************************************************************************/
/***********************************************************************
   Two-Group (X-ray + UV) Flux-Limited Diffusion Implicit Problem Class
   Destructor routine

   author: Daniel R. Reynolds
   date:   January 2013

   PURPOSE: Frees all memory allocated for the DualFLD problem.
************************************************************************/
#ifdef RAD_HYDRO
#include "DualFLD.h"


DualFLD::~DualFLD() {

  // delete HYPRE objects
  HYPRE_StructStencilDestroy(stencil);
  HYPRE_StructGridDestroy(grid);

  // delete EnzoVectors and other internal arrays
  //   EnzoVectors require deleting the structure
  delete sol;
  delete U0;
  delete extsrc;
  delete opacity;

  //   arrays require deleting the array
  HYPRE_StructVectorDestroy(rhsvec);
  HYPRE_StructVectorDestroy(solvec);
  HYPRE_StructMatrixDestroy(P);
  delete[] matentries;
  delete[] rhsentries;
  delete[] HYPREbuff;

  // delete boundary condition arrays
  for (int i=0; i<3; i++)
    for (int j=0; j<2; j++) {
      if (XrBdryVals[i][j] != NULL)
	delete[] XrBdryVals[i][j];
      if (UVBdryVals[i][j] != NULL)
	delete[] UVBdryVals[i][j];
    }
}

#endif
