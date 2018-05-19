/***********************************************************************
/
/  GRID CLASS (WRAPPER FOR EULER SOLVER)
/
/  written by: John Wise
/  date:       May, 2007
/  modified1:  June 1st 2010 E2.0 version
/  modified2:  September 23rd 2010
/              OpenMP loops
/
/  PURPOSE:
/
/  RETURNS:
/    SUCCESS or FAIL
/
************************************************************************/

// Solve the hydro equations with the solver, saving the subgrid fluxes
//


#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
//#include "fortran.def"

int grid::SolvePPM_DE(int CycleNumber, int NumberOfSubgrids, 
		      fluxes *SubgridFluxes[], float *CellWidthTemp[], 
		      Elong_int GridGlobalStart[], int GravityOn, 
		      int NumberOfColours, int colnum[])
{

  int nxz, nyz, nzz, ixyz, nmax;

  nxz = GridEndIndex[0] - GridStartIndex[0] + 1;
  nyz = GridEndIndex[1] - GridStartIndex[1] + 1;
  nzz = GridEndIndex[2] - GridStartIndex[2] + 1;

//  nmax = max(GridDimension[0], GridDimension[1]);
//  nmax = max(GridDimension[2], nmax);

  nmax = 1031;

  ixyz = CycleNumber % GridRank;

  int i,j,k,n;

  int *status = new int[nmax];


  for (n = ixyz; n < ixyz+GridRank; n++) {

    // Update in x-direction
    if ((n % GridRank == 0) && nxz > 1) {

#ifdef INTEL_OMP_SYNTAX
#pragma omp parallel private(k) shared(NumberOfSubgrids, SubgridFluxes, GridGlobalStart, CellWidthTemp, GravityOn, NumberOfColours, colnum, status) default(none)
#else
#pragma omp parallel private(k) shared(GridDimension, NumberOfSubgrids, SubgridFluxes, GridGlobalStart, CellWidthTemp, GravityOn, NumberOfColours, colnum, status) default(none)
#endif
  {

#pragma omp for schedule(static)
      for (k = 0; k < GridDimension[2]; k++) {
       status[k] = 0;
      }

#pragma omp for schedule(dynamic)
      for (k = 0; k < GridDimension[2]; k++) {
        status[k] = this->xEulerSweep(k, NumberOfSubgrids, SubgridFluxes, 
			      GridGlobalStart, CellWidthTemp, GravityOn, 
			      NumberOfColours, colnum);
      } // ENDFOR k

  } // end parallel

    } // ENDIF x-direction

    // Update in y-direction
    if ((n % GridRank == 1) && nyz > 1) {

#ifdef INTEL_OMP_SYNTAX
#pragma omp parallel private(i) shared(NumberOfSubgrids, SubgridFluxes, GridGlobalStart, CellWidthTemp, GravityOn, NumberOfColours, colnum, status) default(none)
#else
#pragma omp parallel private(i) shared(GridDimension, NumberOfSubgrids, SubgridFluxes, GridGlobalStart, CellWidthTemp, GravityOn, NumberOfColours, colnum, status) default(none)
#endif
  {

#pragma omp for schedule(static)
      for (i = 0; i < GridDimension[0]; i++) {
       status[i] = 0;
      }

#pragma omp for schedule(dynamic)
      for (i = 0; i < GridDimension[0]; i++) {
        status[i] = this->yEulerSweep(i, NumberOfSubgrids, SubgridFluxes, 
			      GridGlobalStart, CellWidthTemp, GravityOn, 
			      NumberOfColours, colnum);
      } // ENDFOR i

  } // end parallel

    } // ENDIF y-direction

    // Update in z-direction
    if ((n % GridRank == 2) && nzz > 1) {

#ifdef INTEL_OMP_SYNTAX
#pragma omp parallel private(j) shared(NumberOfSubgrids, SubgridFluxes, GridGlobalStart, CellWidthTemp, GravityOn, NumberOfColours, colnum, status) default(none)
#else
#pragma omp parallel private(j) shared(GridDimension, NumberOfSubgrids, SubgridFluxes, GridGlobalStart, CellWidthTemp, GravityOn, NumberOfColours, colnum, status) default(none)
#endif
  {

#pragma omp for schedule(static)
      for (j = 0; j < GridDimension[1]; j++) {
       status[j] = 0;
      }

#pragma omp for schedule(dynamic)
      for (j = 0; j < GridDimension[1]; j++) {
        status[j] = this->zEulerSweep(j, NumberOfSubgrids, SubgridFluxes, 
			      GridGlobalStart, CellWidthTemp, GravityOn, 
			      NumberOfColours, colnum);
      } // ENDFOR j

  } // end parallel

    } // ENDIF z-direction

  } // ENDFOR n

  delete [] status;

  return SUCCESS;

}
