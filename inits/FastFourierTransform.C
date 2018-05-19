/***********************************************************************
/
/  COMPUTE A FAST FOURIER TRANSFORM (FORWARD OR INVERSE)
/
/  written by: Greg Bryan
/  date:       March, 1995
/  modified1:  Robert Harkness
/  date:       May, 2003
/
/  PURPOSE:
/
/  INPUTS:
/      buffer - field to be FFTed
/      Rank   - rank of FFT
/      DimensionReal[] - declared dimensions of buffer
/      Dimension[]     - active dimensions of buffer
/      direction       - +1 forward, -1 inverse
/
************************************************************************/
 
#include <stdlib.h>
#include <stdio.h>
#include "macros_and_parameters.h"
 
 
int FastFourierTransformPrepareComplex(FLOAT *buffer, int Rank, int DimensionReal[],
                                     int Dimension[], int direction, int type);
 
 
 
 
int FastFourierTransform(FLOAT *buffer, int Rank, int DimensionReal[],
			 int Dimension[], int direction, int type)
{
 

  // A real-to-complex FFT call could be placed here.
  // For simplicity this sets up for a complex-to-complex FFT.
  // FastFourierTransformPrepareComplex is also a shell,
  // the actual choice of a specific FFT routine is 
 
  if (FastFourierTransformPrepareComplex(buffer, Rank, DimensionReal,
				         Dimension, direction, type) == FAIL) {
    fprintf(stderr, "Error in FastFourierTransformPrepareComplex.\n");
    return FAIL;
  }
 
 
  return SUCCESS;
}
