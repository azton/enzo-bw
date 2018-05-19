/***********************************************************************
/
/  CHECK FOR OUTPUT
/
/  written by: Greg Bryan
/  date:       January, 1996
/  modified:   Robert Harkness
/  date:       January, 2007
/              Mods for group in-core i/o
/  date:       May, 2008
/              Remove Dan Reynold's iso_grav code
/
/  PURPOSE:
/    This routine checks a number of criteria for output and then calls
/      the appropriate routine.
/
************************************************************************/

#ifdef RAD_HYDRO
#include "ImplicitProblemABC_preincludes.h"
#endif
 
#include <stdio.h>

#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "Hierarchy.h"
#include "TopGridData.h"
#include "LevelHierarchy.h"
#include "CosmologyParameters.h"
#ifdef RAD_HYDRO
#include "ImplicitProblemABC.h"
#endif
 
/* function prototypes */
 
int WriteAllData(char *basename, int filenumber, HierarchyEntry *TopGrid,
		 TopGridData &MetaData, ExternalBoundary *Exterior,
#ifdef RAD_HYDRO
		 ImplicitProblemABC *ImplicitSolver,
#endif
		 FLOAT WriteTime = -1);

int Group_WriteAllData(char *basename, int filenumber, HierarchyEntry *TopGrid,
		       TopGridData &MetaData, ExternalBoundary *Exterior,
#ifdef RAD_HYDRO
		       ImplicitProblemABC *ImplicitSolver,
#endif
		       FLOAT WriteTime = -1);




int CheckForOutput(HierarchyEntry *TopGrid, TopGridData &MetaData,
		   ExternalBoundary *Exterior,
#ifdef RAD_HYDRO
		   ImplicitProblemABC *ImplicitSolver,
#endif
		   int &WroteData)
{
 
  /* Declarations. */
 
  char *Name;
  int i, Number;
  WroteData = FALSE;
 
  /* Check for output: time-based. */
 
  if (MetaData.Time >= MetaData.TimeLastDataDump + MetaData.dtDataDump
      && MetaData.dtDataDump > 0.0) {
    MetaData.TimeLastDataDump += MetaData.dtDataDump;

#ifdef USE_HDF5_GROUPS
    if (Group_WriteAllData(MetaData.DataDumpName, MetaData.DataDumpNumber++,
		     TopGrid, MetaData, Exterior
#ifdef RAD_HYDRO
		     , ImplicitSolver
#endif
		     ) == FAIL) {
	fprintf(stderr, "Error in Group_WriteAllData.\n");
	return FAIL;
    }
#else
    if (WriteAllData(MetaData.DataDumpName, MetaData.DataDumpNumber++,
		     TopGrid, MetaData, Exterior
#ifdef RAD_HYDRO
		     , ImplicitSolver
#endif
		     ) == FAIL) {
	fprintf(stderr, "Error in WriteAllData.\n");
	return FAIL;
    }
#endif

    WroteData = TRUE;
  }
 
  /* Check for output: cycle-based. */
 
  if (MetaData.CycleNumber >= MetaData.CycleLastDataDump +
                              MetaData.CycleSkipDataDump   &&
      MetaData.CycleSkipDataDump > 0) {
    MetaData.CycleLastDataDump += MetaData.CycleSkipDataDump;

#ifdef USE_HDF5_GROUPS
    if (Group_WriteAllData(MetaData.DataDumpName, MetaData.DataDumpNumber++,
		     TopGrid, MetaData, Exterior
#ifdef RAD_HYDRO
		     , ImplicitSolver
#endif
		     ) == FAIL) {
	fprintf(stderr, "Error in Group_WriteAllData.\n");
	return FAIL;
    }
#else
    if (WriteAllData(MetaData.DataDumpName, MetaData.DataDumpNumber++,
		     TopGrid, MetaData, Exterior
#ifdef RAD_HYDRO
		     , ImplicitSolver
#endif               
		     ) == FAIL) {
	fprintf(stderr, "Error in WriteAllData.\n");
	return FAIL;
    }
#endif

    WroteData = TRUE;
  }
 
  /* Check for output: redshift-based. */
 
  if (ComovingCoordinates)
    for (i = 0; i < MAX_NUMBER_OF_OUTPUT_REDSHIFTS; i++)
      if (CosmologyOutputRedshift[i] != -1)
	if (MetaData.Time >= CosmologyOutputRedshiftTime[i]) {
	  CosmologyOutputRedshift[i] = -1; // done, turn it off
	  if (CosmologyOutputRedshiftName[i] == NULL) {
	    Name   = MetaData.RedshiftDumpName;
	    Number = i;   // append number to end of name
	  }
	  else {
	    Name   = CosmologyOutputRedshiftName[i];
	    Number = -1;  // Don't append number (####) to end of name
	  }

#ifdef USE_HDF5_GROUPS
	  if (Group_WriteAllData(Name, Number, TopGrid, MetaData, Exterior
#ifdef RAD_HYDRO
				 , ImplicitSolver
#endif
				 ) == FAIL) {
	    fprintf(stderr, "Error in Group_WriteAllData.\n");
	    return FAIL;
	  }
#else
	  if (WriteAllData(Name, Number, TopGrid, MetaData, Exterior
#ifdef RAD_HYDRO
			   , ImplicitSolver
#endif
			   ) == FAIL) {
	    fprintf(stderr, "Error in WriteAllData.\n");
	    return FAIL;
	  }
#endif

	  WroteData = TRUE;
	}
 
  return SUCCESS;
}
