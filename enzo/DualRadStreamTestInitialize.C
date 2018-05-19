/***********************************************************************
/
/  Initialization routine for streaming UV + X-ray radiation problem
/
/  written by: Daniel Reynolds
/  date:       April 2013
/
/  PURPOSE:
/
/  RETURNS: SUCCESS or FAIL
/
************************************************************************/

// This routine intializes a new simulation based on the parameter file.
#include <string.h>
#include <stdio.h>
#include <math.h>
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "Hierarchy.h"
#include "TopGridData.h"


/* default constants */
#define DEFAULT_MU 0.6       // mean molecular mass
#define MIN_TEMP 1.0         // minimum temperature [K]


// function prototypes
int InitializeRateData(FLOAT Time);



int DualRadStreamTestInitialize(FILE *fptr, FILE *Outfptr, 
				HierarchyEntry &TopGrid,
				TopGridData &MetaData, int local) 
{

#ifdef RAD_HYDRO

  // if (debug)
  //   fprintf(stdout,"Entering DualRadStreamTestInitialize routine\n");

  char *kphHIName    = "HI_kph";
  char *kphHeIName   = "HeI_kph";
  char *kphHeIIName  = "HeII_kph";
  char *gammaName    = "PhotoGamma";
  char *kdissH2IName = "H2I_kdiss";
  char *DensName  = "Density";
  char *TEName    = "Total_Energy";
  char *IEName    = "Internal_Energy";
  char *Vel0Name  = "x-velocity";
  char *Vel1Name  = "y-velocity";
  char *Vel2Name  = "z-velocity";
  char *RadName0  = "Xray_Radiation";
  char *RadName1  = "UV_Radiation";
  char *HIName    = "HI_Density";
  char *HIIName   = "HII_Density";
  char *HeIName   = "HeI_Density";
  char *HeIIName  = "HeII_Density";
  char *HeIIIName = "HeIII_Density";
  char *DeName    = "Electron_Density";

  // local declarations
  char line[MAX_LINE_LENGTH];
  int  dim, ret;

  // Setup and parameters:
  //  1. ambient density (should be very small) - free parameter
  //  2. ambient gas velocity - free parameter
  //  3. ambient gas temperature
  //  4. ambient radiation energy
  //  5. Hydrogen mass fraction 
  //  6. initial fraction HII
  //  7. initial fraction HeII
  //  8. initial fraction HeIII
  //  9. Number of chemical species
  // 10. dimension for streaming test problem {0,1,2}
  // 11. direction for streaming radiation {0:l->r, 1:r->l}
  // 12. type of streaming radiation {0:Xray, 1:UV}
  float X0Velocity           = 0.0;
  float X1Velocity           = 0.0;
  float X2Velocity           = 0.0;
  float Density              = 1.0;
  float Temperature          = 1.0;
  float IEnergy              = -1.0;
  float UVRadiation          = 1.0;
  float XrRadiation          = 1.0;
  float HydrogenMassFraction = 1.0;
  float InitialFractionHII   = 0.0;
  float InitialFractionHeII  = 0.0;
  float InitialFractionHeIII = 0.0;
  int   Chemistry            = 1;
  int   RadStreamDim         = 0;
  int   RadStreamDir         = 0;
  int   RadStreamType        = 0;

  // overwrite parameters from RadHydroParamFile file, if it exists
   if (MetaData.RadHydroParameterFname != NULL) {
    FILE *RHfptr;
    if ((RHfptr = fopen(MetaData.RadHydroParameterFname, "r")) != NULL) {
      while (fgets(line, MAX_LINE_LENGTH, RHfptr) != NULL) {
	ret = 0;
	// read relevant problem parameters
	ret += sscanf(line, "DualFLDChemistry = %"ISYM, &Chemistry);
	ret += sscanf(line, "DualFLDHFraction = %"FSYM, &HydrogenMassFraction);
	ret += sscanf(line, "Velocity = %"FSYM" %"FSYM" %"FSYM,
		      &X0Velocity, &X1Velocity, &X2Velocity);
	ret += sscanf(line, "Density = %"FSYM, &Density);
	ret += sscanf(line, "Temperature = %"FSYM, &Temperature);
	ret += sscanf(line, "IEnergy = %"FSYM, &IEnergy);
	ret += sscanf(line, "UVRadiation = %"FSYM, &UVRadiation);
	ret += sscanf(line, "XrayRadiation = %"FSYM, &XrRadiation);
	ret += sscanf(line, "InitialFractionHII = %"FSYM, &InitialFractionHII);
	ret += sscanf(line, "InitialFractionHeII = %"FSYM, 
		      &InitialFractionHeII);
	ret += sscanf(line, "InitialFractionHeIII = %"FSYM, 
		      &InitialFractionHeIII);
	ret += sscanf(line, "RadStreamDim = %"ISYM, &RadStreamDim);
	ret += sscanf(line, "RadStreamDir = %"ISYM, &RadStreamDir);
	ret += sscanf(line, "RadStreamType = %"ISYM, &RadStreamType);
      } // end input from parameter file
      fclose(RHfptr);
    }
  }

  // ensure that streaming dimension is active for this rank
  if (RadStreamDim >= MetaData.TopGridRank) {
    fprintf(stderr,"RadStreamDim = %"ISYM" illegal for rank = %"ISYM"!\n",
	    RadStreamDim,MetaData.TopGridRank);
    return FAIL;
  }

  // ensure that streaming direction is allowable
  if (RadStreamDir < 0 || RadStreamDir > 1) {
    fprintf(stderr,"RadStreamDir = %"ISYM" illegal!\n", RadStreamDir);
    return FAIL;
  }

  // ensure that streaming radiation type is allowable
  if (RadStreamType < 0 || RadStreamType > 1) {
    fprintf(stderr,"RadStreamType = %"ISYM" illegal!\n", RadStreamType);
    return FAIL;
  }

  // set up CoolData object if not already set up
  if (CoolData.ceHI == NULL) 
    if (InitializeRateData(MetaData.Time) == FAIL) {
      fprintf(stderr,"Error in InitializeRateData.\n");
      return FAIL;
    }


  // if temperature specified and not internal energy, perform conversion here
  if (IEnergy == -1.0) {
    if (Temperature == -1.0) {
      fprintf(stderr,"Initialize error: either Temperature or IEnergy required!\n");
      return FAIL;
    }
    else {
      Temperature = max(Temperature, MIN_TEMP);  // enforce minimum
      float mp = 1.67262171e-24;    // proton mass [g]
      float kb = 1.3806504e-16;     // boltzmann constant [erg/K]
      float nH, HI, HII, nHe, HeI, HeII, HeIII, ne, num_dens, mu;
      if (Chemistry == 1) {
	nH = Density*HydrogenMassFraction;
	HI = nH*(1.0 - InitialFractionHII);
	HII = nH*InitialFractionHII;
	ne = HII;
	num_dens = HI + HII + ne;
	mu = Density/num_dens;
      }
      else if (Chemistry == 3) {
	nH = Density*HydrogenMassFraction;
	nHe = Density*(1.0 - HydrogenMassFraction);
	HI = nH*(1.0 - InitialFractionHII);
	HII = nH*InitialFractionHII;
	HeII = nHe*InitialFractionHeII;
	HeIII = nHe*InitialFractionHeIII;
	HeI = nHe - HeII - HeIII;
	ne = HII + HeII/4.0 + HeIII/2.0;
	num_dens = 0.25*(HeI + HeII + HeIII) + HI + HII + ne;
	mu = Density/num_dens;
      }
      else {
	fprintf(stderr,"Initialize error: DualFLDChemistry != {1,3}\n");
	return FAIL;	
      }
      // compute the internal energy
      IEnergy = kb*Temperature/mu/mp/(Gamma-1.0);	
    }
  }

  // set up the grid(s) on this level
  HierarchyEntry *Temp = &TopGrid;
  while (Temp != NULL) {
    if (Temp->GridData->DualRadStreamTestInitializeGrid(
		        Chemistry, Density, X0Velocity, X1Velocity, 
			X2Velocity, IEnergy, UVRadiation, XrRadiation, 
			HydrogenMassFraction, InitialFractionHII, 
			InitialFractionHeII, InitialFractionHeIII, 
			RadStreamDim, RadStreamDir, RadStreamType, 
			local) == FAIL) {
      fprintf(stderr, "Error in DualRadStreamTestInitializeGrid.\n");
      return FAIL;
    }
    Temp = Temp->NextGridThisLevel;
  }

  // set up field names and units
  int BaryonField = 0;
  DataLabel[BaryonField++] = DensName;
  DataLabel[BaryonField++] = TEName;
  if (DualEnergyFormalism) 
    DataLabel[BaryonField++] = IEName;
  DataLabel[BaryonField++] = Vel0Name;
  DataLabel[BaryonField++] = Vel1Name;
  DataLabel[BaryonField++] = Vel2Name;
  DataLabel[BaryonField++] = RadName0;
  DataLabel[BaryonField++] = RadName1;
  DataLabel[BaryonField++] = DeName;
  DataLabel[BaryonField++] = HIName;
  DataLabel[BaryonField++] = HIIName;
  if ((Chemistry == 3) || (MultiSpecies > 0)) {
    DataLabel[BaryonField++] = HeIName;
    DataLabel[BaryonField++] = HeIIName;
    DataLabel[BaryonField++] = HeIIIName;
  }

  // if using external chemistry/cooling, set rate labels
  if (RadiativeCooling) {
    DataLabel[BaryonField++] = kphHIName;
    DataLabel[BaryonField++] = gammaName;
    if (RadiativeTransferHydrogenOnly == FALSE) {
      DataLabel[BaryonField++] = kphHeIName;
      DataLabel[BaryonField++] = kphHeIIName;
    }
    if (MultiSpecies > 1)
      DataLabel[BaryonField++] = kdissH2IName;
  }

  for (int i=0; i<BaryonField; i++) 
    DataUnits[i] = NULL;

  return SUCCESS;

#else

  fprintf(stderr,"Error: RAD_HYDRO must be enabled for this test!\n");
  return FAIL;

#endif

}

