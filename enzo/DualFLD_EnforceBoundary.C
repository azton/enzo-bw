/*****************************************************************************
 * Copyright 2013 Daniel R. Reynolds                                         *
 *                                                                           *
 * This software is released under the terms of the "Enzo Public License"    *
 * in the accompanying LICENSE file.                                         *
 *****************************************************************************/
/***********************************************************************
   Two-Group (X-ray + UV) Flux-Limited Diffusion Implicit Problem Class
   EnforceBoundary routine

   author: Daniel R. Reynolds
   date:   January 2013

   PURPOSE: Enforces boundary conditions on a FLD problem vector.

      Note: Neumann values are enforced on the first layer of ghost
            zones using a first-order central difference approximation 
            to the first (outward-normal) derivative.
      Note: Since the internal radiation variables are comoving and
            normalized, we renormalize the boundary conditions as they 
            are enforced to match the internal units.
************************************************************************/
#ifdef RAD_HYDRO
#include "DualFLD.h"


int DualFLD::EnforceBoundary(EnzoVector *u) {

  // get information about the vector u, and check against BC dims
  int i, i2, j, j2, k, k2, idx, idx2, idxbc, ifield;
  int udims[4], ugh[3][2];
  u->size(&udims[0], &udims[1], &udims[2], &udims[3], 
	  &ugh[0][0], &ugh[0][1], &ugh[1][0], 
	  &ugh[1][1], &ugh[2][0], &ugh[2][1]);
  if (udims[0] != LocDims[0]) {
    fprintf(stderr,"p%"ISYM" EnforceBC: mismatched x0 dims %"ISYM"!=%"ISYM"\n",
	    MyProcessorNumber,udims[0],LocDims[0]);
    return FAIL;
  }
  if (udims[1] != LocDims[1]) {
    fprintf(stderr,"p%"ISYM" EnforceBC: mismatched x1 dims %"ISYM"!=%"ISYM"\n",
	    MyProcessorNumber,udims[1],LocDims[1]);
    return FAIL;
  }
  if (udims[2] != LocDims[2]) {
    fprintf(stderr,"p%"ISYM" EnforceBC: mismatched x2 dims %"ISYM"!=%"ISYM"\n",
	    MyProcessorNumber,udims[2],LocDims[2]);
    return FAIL;
  }

  // set some shortcuts, data accessor pointer
  float *udata;
  int x0len = udims[0] + ugh[0][0] + ugh[0][1];
  int x1len = udims[1] + ugh[1][0] + ugh[1][1];
  float dxa, dya, dza, EUnits;
  dxa = dx[0]*LenUnits;
  if (rank > 1)  dya = dx[1]*LenUnits;
  if (rank > 2)  dza = dx[2]*LenUnits;

  
  // enforce boundary conditions on Xray radiation field
  udata = u->GetData(0);
  EUnits = XrUnits;

  // x0 left boundary
  //   Dirichlet
  if (OnBdry[0][0] && (XrBdryType[0][0]==1)) {
    for (k=0; k<LocDims[2]; k++)
      for (j=0; j<LocDims[1]; j++)
	for (i=0; i<ugh[0][0]; i++) {
	  idxbc = k*LocDims[1] + j;
	  idx = ((k+ugh[2][0])*x1len + j+ugh[1][0])*x0len + i;
	  udata[idx] = XrBdryVals[0][0][idxbc]/EUnits;
	}
  }
  //   Neumann
  if (OnBdry[0][0] && (XrBdryType[0][0]==2)) {
    i = -1;  i2 = i+1;
    for (k=0; k<LocDims[2]; k++)
      for (j=0; j<LocDims[1]; j++) {
	idx = ((k+ugh[2][0])*x1len + j+ugh[1][0])*x0len + i+ugh[0][0];
	idx2 = ((k+ugh[2][0])*x1len + j+ugh[1][0])*x0len + i2+ugh[0][0];
	idxbc = k*LocDims[1] + j;
	udata[idx] = udata[idx2] + dxa*XrBdryVals[0][0][idxbc]/EUnits;
      }
  }

  // x0 right boundary
  //   Dirichlet
  if (OnBdry[0][1] && (XrBdryType[0][1]==1)) {
    for (k=0; k<LocDims[2]; k++)
      for (j=0; j<LocDims[1]; j++)
	for (i=ArrDims[0]-ugh[0][1]; i<ArrDims[0]; i++) {
	  idxbc = k*LocDims[1] + j;
	  idx = ((k+ugh[2][0])*x1len + j+ugh[1][0])*x0len + i;
	  udata[idx] = XrBdryVals[0][1][idxbc]/EUnits;
	}
  }
  //   Neumann
  if (OnBdry[0][1] && (XrBdryType[0][1]==2)) {
    i = LocDims[0];  i2 = i-1;
    for (k=0; k<LocDims[2]; k++)
      for (j=0; j<LocDims[1]; j++) {
	idx = ((k+ugh[2][0])*x1len + j+ugh[1][0])*x0len + i+ugh[0][0];
	idx2 = ((k+ugh[2][0])*x1len + j+ugh[1][0])*x0len + i2+ugh[0][0];
	idxbc = k*LocDims[1] + j;
	udata[idx] = udata[idx2] + dxa*XrBdryVals[0][1][idxbc]/EUnits;
      }
  }

  if (rank > 1) {
    // x1 left boundary
    //   Dirichlet
    if (OnBdry[1][0] && (XrBdryType[1][0]==1)) {
      for (k=0; k<LocDims[2]; k++)
	for (j=0; j<ugh[1][0]; j++)
	  for (i=0; i<LocDims[0]; i++) {
	    idx = ((k+ugh[2][0])*x1len + j)*x0len + i+ugh[0][0];
	    idxbc = i*LocDims[2] + k;
	    udata[idx] = XrBdryVals[1][0][idxbc]/EUnits;
	  }
    }
    //   Neumann
    if (OnBdry[1][0] && (XrBdryType[1][0]==2)) {
      j = -1;  j2 = j+1;
      for (k=0; k<LocDims[2]; k++)
	for (i=0; i<LocDims[0]; i++) {
	  idx = ((k+ugh[2][0])*x1len + j+ugh[1][0])*x0len + i+ugh[0][0];
	  idx2 = ((k+ugh[2][0])*x1len + j2+ugh[1][0])*x0len + i+ugh[0][0];
	  idxbc = i*LocDims[2] + k;
	  udata[idx] = udata[idx2] + dya*XrBdryVals[1][0][idxbc]/EUnits;
	}
    }
      
    // x1 right boundary
    //   Dirichlet
    if (OnBdry[1][1] && (XrBdryType[1][1]==1)) {
      for (k=0; k<LocDims[2]; k++)
	for (j=ArrDims[1]-ugh[1][1]; j<ArrDims[1]; j++)
	  for (i=0; i<LocDims[0]; i++) {
	    idx = ((k+ugh[2][0])*x1len + j)*x0len + i+ugh[0][0];
	    idxbc = i*LocDims[2] + k;
	    udata[idx] = XrBdryVals[1][1][idxbc]/EUnits;
	  }
    }
    //   Neumann
    if (OnBdry[1][1] && (XrBdryType[1][1]==2)) {
      j = LocDims[1];  j2 = j-1;
      for (k=0; k<LocDims[2]; k++)
	for (i=0; i<LocDims[0]; i++) {
	  idx = ((k+ugh[2][0])*x1len + j+ugh[1][0])*x0len + i+ugh[0][0];
	  idx2 = ((k+ugh[2][0])*x1len + j2+ugh[1][0])*x0len + i+ugh[0][0];
	  idxbc = i*LocDims[2] + k;
	  udata[idx] = udata[idx2] + dya*XrBdryVals[1][1][idxbc]/EUnits;
	}
    }
  }  // end if rank > 1
     
  if (rank > 2) {
    // x2 left boundary
    //   Dirichlet
    if (OnBdry[2][0] && (XrBdryType[2][0]==1)) {
      for (k=0; k<ugh[2][0]; k++)
	for (j=0; j<LocDims[1]; j++)
	  for (i=0; i<LocDims[0]; i++) {
	    idx = (k*x1len + j+ugh[1][0])*x0len + i+ugh[0][0];
	    idxbc = j*LocDims[0] + i;
	    udata[idx] = XrBdryVals[2][0][idxbc]/EUnits;
	  }
    }
    //   Neumann
    if (OnBdry[2][0] && (XrBdryType[2][0]==2)) {
      k = -1;  k2 = k+1;
      for (j=0; j<LocDims[1]; j++)
	for (i=0; i<LocDims[0]; i++) {
	  idx = ((k+ugh[2][0])*x1len + j+ugh[1][0])*x0len + i+ugh[0][0];
	  idx2 = ((k2+ugh[2][0])*x1len + j+ugh[1][0])*x0len + i+ugh[0][0];
	  idxbc = j*LocDims[0] + i;
	  udata[idx] = udata[idx2] + dza*XrBdryVals[2][0][idxbc]/EUnits;
	}
    }
      
    // x2 right boundary
    //   Dirichlet
    if (OnBdry[2][1] && (XrBdryType[2][1]==1)) {
      for (k=ArrDims[2]-ugh[2][1]; k<ArrDims[2]; k++)
	for (j=0; j<LocDims[1]; j++)
	  for (i=0; i<LocDims[0]; i++) {
	    idx = (k*x1len + j+ugh[1][0])*x0len + i+ugh[0][0];
	    idxbc = j*LocDims[0] + i;
	    udata[idx] = XrBdryVals[2][1][idxbc]/EUnits;
	  }
    }
    //   Neumann
    if (OnBdry[2][1] && (XrBdryType[2][1]==2)) {
      k = LocDims[2];  k2 = k-1;
      for (j=0; j<LocDims[1]; j++)
	for (i=0; i<LocDims[0]; i++) {
	  idx = ((k+ugh[2][0])*x1len + j+ugh[1][0])*x0len + i+ugh[0][0];
	  idx2 = ((k2+ugh[2][0])*x1len + j+ugh[1][0])*x0len + i+ugh[0][0];
	  idxbc = j*LocDims[0] + i;
	  udata[idx] = udata[idx2] + dza*XrBdryVals[2][1][idxbc]/EUnits;
	}
    }
  }  // end if rank > 2



  if (!XrayOnly) {

    // enforce boundary conditions on UV radiation field
    udata = u->GetData(1);
    EUnits = UVUnits;

    // x0 left boundary
    //   Dirichlet
    if (OnBdry[0][0] && (UVBdryType[0][0]==1)) {
      for (k=0; k<LocDims[2]; k++)
	for (j=0; j<LocDims[1]; j++)
	  for (i=0; i<ugh[0][0]; i++) {
	    idxbc = k*LocDims[1] + j;
	    idx = ((k+ugh[2][0])*x1len + j+ugh[1][0])*x0len + i;
	    udata[idx] = UVBdryVals[0][0][idxbc]/EUnits;
	  }
    }
    //   Neumann
    if (OnBdry[0][0] && (UVBdryType[0][0]==2)) {
      i = -1;  i2 = i+1;
      for (k=0; k<LocDims[2]; k++)
	for (j=0; j<LocDims[1]; j++) {
	  idx = ((k+ugh[2][0])*x1len + j+ugh[1][0])*x0len + i+ugh[0][0];
	  idx2 = ((k+ugh[2][0])*x1len + j+ugh[1][0])*x0len + i2+ugh[0][0];
	  idxbc = k*LocDims[1] + j;
	  udata[idx] = udata[idx2] + dxa*UVBdryVals[0][0][idxbc]/EUnits;
	}
    }

    // x0 right boundary
    //   Dirichlet
    if (OnBdry[0][1] && (UVBdryType[0][1]==1)) {
      for (k=0; k<LocDims[2]; k++)
	for (j=0; j<LocDims[1]; j++)
	  for (i=ArrDims[0]-ugh[0][1]; i<ArrDims[0]; i++) {
	    idxbc = k*LocDims[1] + j;
	    idx = ((k+ugh[2][0])*x1len + j+ugh[1][0])*x0len + i;
	    udata[idx] = UVBdryVals[0][1][idxbc]/EUnits;
	  }
    }
    //   Neumann
    if (OnBdry[0][1] && (UVBdryType[0][1]==2)) {
      i = LocDims[0];  i2 = i-1;
      for (k=0; k<LocDims[2]; k++)
	for (j=0; j<LocDims[1]; j++) {
	  idx = ((k+ugh[2][0])*x1len + j+ugh[1][0])*x0len + i+ugh[0][0];
	  idx2 = ((k+ugh[2][0])*x1len + j+ugh[1][0])*x0len + i2+ugh[0][0];
	  idxbc = k*LocDims[1] + j;
	  udata[idx] = udata[idx2] + dxa*UVBdryVals[0][1][idxbc]/EUnits;
	}
    }

    if (rank > 1) {
      // x1 left boundary
      //   Dirichlet
      if (OnBdry[1][0] && (UVBdryType[1][0]==1)) {
	for (k=0; k<LocDims[2]; k++)
	  for (j=0; j<ugh[1][0]; j++)
	    for (i=0; i<LocDims[0]; i++) {
	      idx = ((k+ugh[2][0])*x1len + j)*x0len + i+ugh[0][0];
	      idxbc = i*LocDims[2] + k;
	      udata[idx] = UVBdryVals[1][0][idxbc]/EUnits;
	    }
      }
      //   Neumann
      if (OnBdry[1][0] && (UVBdryType[1][0]==2)) {
	j = -1;  j2 = j+1;
	for (k=0; k<LocDims[2]; k++)
	  for (i=0; i<LocDims[0]; i++) {
	    idx = ((k+ugh[2][0])*x1len + j+ugh[1][0])*x0len + i+ugh[0][0];
	    idx2 = ((k+ugh[2][0])*x1len + j2+ugh[1][0])*x0len + i+ugh[0][0];
	    idxbc = i*LocDims[2] + k;
	    udata[idx] = udata[idx2] + dya*UVBdryVals[1][0][idxbc]/EUnits;
	  }
      }
      
      // x1 right boundary
      //   Dirichlet
      if (OnBdry[1][1] && (UVBdryType[1][1]==1)) {
	for (k=0; k<LocDims[2]; k++)
	  for (j=ArrDims[1]-ugh[1][1]; j<ArrDims[1]; j++)
	    for (i=0; i<LocDims[0]; i++) {
	      idx = ((k+ugh[2][0])*x1len + j)*x0len + i+ugh[0][0];
	      idxbc = i*LocDims[2] + k;
	      udata[idx] = UVBdryVals[1][1][idxbc]/EUnits;
	    }
      }
      //   Neumann
      if (OnBdry[1][1] && (UVBdryType[1][1]==2)) {
	j = LocDims[1];  j2 = j-1;
	for (k=0; k<LocDims[2]; k++)
	  for (i=0; i<LocDims[0]; i++) {
	    idx = ((k+ugh[2][0])*x1len + j+ugh[1][0])*x0len + i+ugh[0][0];
	    idx2 = ((k+ugh[2][0])*x1len + j2+ugh[1][0])*x0len + i+ugh[0][0];
	    idxbc = i*LocDims[2] + k;
	    udata[idx] = udata[idx2] + dya*UVBdryVals[1][1][idxbc]/EUnits;
	  }
      }
    }  // end if rank > 1
     
    if (rank > 2) {
      // x2 left boundary
      //   Dirichlet
      if (OnBdry[2][0] && (UVBdryType[2][0]==1)) {
	for (k=0; k<ugh[2][0]; k++)
	  for (j=0; j<LocDims[1]; j++)
	    for (i=0; i<LocDims[0]; i++) {
	      idx = (k*x1len + j+ugh[1][0])*x0len + i+ugh[0][0];
	      idxbc = j*LocDims[0] + i;
	      udata[idx] = UVBdryVals[2][0][idxbc]/EUnits;
	    }
      }
      //   Neumann
      if (OnBdry[2][0] && (UVBdryType[2][0]==2)) {
	k = -1;  k2 = k+1;
	for (j=0; j<LocDims[1]; j++)
	  for (i=0; i<LocDims[0]; i++) {
	    idx = ((k+ugh[2][0])*x1len + j+ugh[1][0])*x0len + i+ugh[0][0];
	    idx2 = ((k2+ugh[2][0])*x1len + j+ugh[1][0])*x0len + i+ugh[0][0];
	    idxbc = j*LocDims[0] + i;
	    udata[idx] = udata[idx2] + dza*UVBdryVals[2][0][idxbc]/EUnits;
	  }
      }
      
      // x2 right boundary
      //   Dirichlet
      if (OnBdry[2][1] && (UVBdryType[2][1]==1)) {
	for (k=ArrDims[2]-ugh[2][1]; k<ArrDims[2]; k++)
	  for (j=0; j<LocDims[1]; j++)
	    for (i=0; i<LocDims[0]; i++) {
	      idx = (k*x1len + j+ugh[1][0])*x0len + i+ugh[0][0];
	      idxbc = j*LocDims[0] + i;
	      udata[idx] = UVBdryVals[2][1][idxbc]/EUnits;
	    }
      }
      //   Neumann
      if (OnBdry[2][1] && (UVBdryType[2][1]==2)) {
	k = LocDims[2];  k2 = k-1;
	for (j=0; j<LocDims[1]; j++)
	  for (i=0; i<LocDims[0]; i++) {
	    idx = ((k+ugh[2][0])*x1len + j+ugh[1][0])*x0len + i+ugh[0][0];
	    idx2 = ((k2+ugh[2][0])*x1len + j+ugh[1][0])*x0len + i+ugh[0][0];
	    idxbc = j*LocDims[0] + i;
	    udata[idx] = udata[idx2] + dza*UVBdryVals[2][1][idxbc]/EUnits;
	  }
      }
    }  // end if rank > 2

  }  // end if (!XrayOnly)
      
  return SUCCESS;
}

#endif
