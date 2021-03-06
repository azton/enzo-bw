/***********************************************************************
/
/  GRID CLASS (INITIALIZE THE RAD-HYDRO CLUMP IONIZATION TEST) 
/
/  written by: Daniel Reynolds
/  date:       December 2007
/
/  PURPOSE: 
/
/  RETURNS: FAIL or SUCCESS
/
************************************************************************/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"


// function prototypes
int GetUnits(float *DensityUnits, float *LengthUnits, 
	     float *TemperatureUnits, float *TimeUnits, 
	     float *VelocityUnits, float *MassUnits, FLOAT Time);



int grid::RHIonizationClumpInitializeGrid(int NumChemicals, 
					  float NumDensityIn, 
					  float NumDensityOut, 
					  float VxConst, 
					  float VyConst, 
					  float VzConst, 
					  float IEConstIn, 
					  float IEConstOut, 
					  float EgConst, 
					  float HMassFrac, 
					  float InitFracHII, 
					  float InitFracHeII, 
					  float InitFracHeIII, 
					  float ClumpCenterX,
					  float ClumpCenterY,
					  float ClumpCenterZ,
					  float ClumpRadius,
					  int   local)
{
#ifdef RAD_HYDRO

//   if (debug)
//     fprintf(stdout,"Entering grid::RHIonizationClumpInitializeGrid routine\n");

  // determine whether data should be allocated/initialized or not
  int NewData = TRUE;
  if ((ParallelRootGridIO == TRUE) && (local == 0))
    NewData = FALSE;

  // if grids allocated and already set up (i.e. restart), return
  if ((NumberOfBaryonFields > 5) && (BaryonField[5] != NULL))
    return SUCCESS;


  // create necessary baryon fields
  int RhoNum, TENum, IENum, V0Num, V1Num, V2Num, EgNum, DeNum, 
    HINum, HIINum, HeINum, HeIINum, HeIIINum, kphHINum, kphHeINum, 
    kphHeIINum, gammaNum, kdissH2INum;
  NumberOfBaryonFields = 0;
  FieldType[RhoNum = NumberOfBaryonFields++] = Density;
  FieldType[TENum = NumberOfBaryonFields++]  = TotalEnergy;
  if (DualEnergyFormalism) 
    FieldType[IENum = NumberOfBaryonFields++]  = InternalEnergy;
  FieldType[V0Num = NumberOfBaryonFields++]    = Velocity1;
  FieldType[V1Num = NumberOfBaryonFields++]    = Velocity2;
  FieldType[V2Num = NumberOfBaryonFields++]    = Velocity3;
  FieldType[EgNum = NumberOfBaryonFields++]    = RadiationFreq0;
  FieldType[DeNum = NumberOfBaryonFields++]    = ElectronDensity;
  FieldType[HINum = NumberOfBaryonFields++]    = HIDensity;
  FieldType[HIINum = NumberOfBaryonFields++]   = HIIDensity;
  if ((NumChemicals == 3) || (MultiSpecies > 0)) {
    FieldType[HeINum   = NumberOfBaryonFields++] = HeIDensity;
    FieldType[HeIINum  = NumberOfBaryonFields++] = HeIIDensity;
    FieldType[HeIIINum = NumberOfBaryonFields++] = HeIIIDensity;
  }
  // if using external chemistry/cooling, set rate fields
  if (RadiativeCooling) {
    FieldType[kphHINum = NumberOfBaryonFields++] = kphHI;
    FieldType[gammaNum = NumberOfBaryonFields++] = PhotoGamma;
    if (RadiativeTransferHydrogenOnly == FALSE) {
      FieldType[kphHeINum  = NumberOfBaryonFields++] = kphHeI;
      FieldType[kphHeIINum = NumberOfBaryonFields++] = kphHeII;
    }
    if (MultiSpecies > 1)
      FieldType[kdissH2INum = NumberOfBaryonFields++] = kdissH2I;
  }


  // set the subgrid static flag (necessary??)
  SubgridsAreStatic = FALSE;

  // Return if this doesn't concern us.
  if (ProcessorNumber != MyProcessorNumber)
    return SUCCESS;

  // Get various units
  float DensityUnits=1, LengthUnits=1, TemperatureUnits=1, TimeUnits=1,
    VelocityUnits=1, MassUnits=1;
  if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	       &TimeUnits, &VelocityUnits, &MassUnits, Time) == FAIL) {
    fprintf(stderr,"Error in GetUnits.\n");
    return FAIL;
  }
  if (debug) {
    fprintf(stdout,"  Internal Unit Conversion Factors:\n");
    fprintf(stdout,"         length = %g\n",LengthUnits);
    fprintf(stdout,"           mass = %g\n",MassUnits);
    fprintf(stdout,"           time = %g\n",TimeUnits);
  }

  // compute size of fields
  int size = 1;
  for (int dim=0; dim<GridRank; dim++)  size *= GridDimension[dim];
 
  // allocate fields
  if (NewData == TRUE) {
    for (int field=0; field<NumberOfBaryonFields; field++)
      if (BaryonField[field] == NULL)
	BaryonField[field] = new float[size];
 
    // set fluid density, total energy, [internal energy,] velocities, 
    // radiation energy, electron density, chemical species
    int i, j, k;
    float mp = 1.67262171e-24;    // proton mass [g]
    float TEConstIn = (IEConstIn + 0.5*(VxConst*VxConst + 
					VyConst*VyConst + 
					VzConst*VzConst));
    float TEConstOut = (IEConstOut + 0.5*(VxConst*VxConst + 
					  VyConst*VyConst + 
					  VzConst*VzConst));
    float RhoConstIn    = NumDensityIn * mp;
    float RhoConstOut   = NumDensityOut * mp;
    float HIIConstIn    = InitFracHII * HMassFrac * RhoConstIn;
    float HIIConstOut   = InitFracHII * HMassFrac * RhoConstOut;
    float HIConstIn     = HMassFrac * RhoConstIn - HIIConstIn;
    float HIConstOut    = HMassFrac * RhoConstOut - HIIConstOut;
    float HeIIConstIn   = InitFracHeII * (1.0 - HMassFrac) * RhoConstIn;
    float HeIIConstOut  = InitFracHeII * (1.0 - HMassFrac) * RhoConstOut;
    float HeIIIConstIn  = InitFracHeIII * (1.0 - HMassFrac) * RhoConstIn;
    float HeIIIConstOut = InitFracHeIII * (1.0 - HMassFrac) * RhoConstOut;
    float HeIConstIn    = (1.0 - HMassFrac) * RhoConstIn - HeIIConstIn - HeIIIConstIn;
    float HeIConstOut   = (1.0 - HMassFrac) * RhoConstOut - HeIIConstOut - HeIIIConstOut;
    float DeConstIn     = HIIConstIn + 0.25*HeIIConstIn + 0.5*HeIIIConstIn;
    float DeConstOut    = HIIConstOut + 0.25*HeIIConstOut + 0.5*HeIIIConstOut;
    float eUnits = VelocityUnits*VelocityUnits;
    float EUnits = DensityUnits*eUnits;
    // initialize clump-independent quantities
    for (i=0; i<size; i++) {
      BaryonField[V0Num][i]  = VxConst/VelocityUnits;
      BaryonField[V1Num][i]  = VyConst/VelocityUnits;
      BaryonField[V2Num][i]  = VzConst/VelocityUnits;
      BaryonField[EgNum][i]  = EgConst/EUnits;
    }
    // if using external chemistry/cooling, set rate fields
    if (RadiativeCooling) {
      for (i=0; i<size; i++)  BaryonField[kphHINum][i] = 0.0;
      for (i=0; i<size; i++)  BaryonField[gammaNum][i] = 0.0;
      if (RadiativeTransferHydrogenOnly == FALSE) {
	for (i=0; i<size; i++)  BaryonField[kphHeINum][i]  = 0.0;
	for (i=0; i<size; i++)  BaryonField[kphHeIINum][i] = 0.0;
      }
      if (MultiSpecies > 1)
	for (i=0; i<size; i++)  BaryonField[kdissH2INum][i] = 0.0;
    }
    
    // initialize clump-dependent quantities
    float Vint, Vext;
    float gridx0l = GridLeftEdge[0];
    float gridx0r = GridRightEdge[0];
    float gridx1l = GridLeftEdge[1];
    float gridx1r = GridRightEdge[1];
    float gridx2l = GridLeftEdge[2];
    float gridx2r = GridRightEdge[2];
    float dx0 = (gridx0r-gridx0l)/GridDimension[0];
    float dx1 = (gridx1r-gridx1l)/GridDimension[1];
    float dx2 = (gridx2r-gridx2l)/GridDimension[2];
    int idx, l;
    float x0l, x0r, x1l, x1r, x2l, x2r, Vin, Vout;
    float d[8];
    for (k=0; k<GridDimension[2]; k++) {
      x2l = gridx2l + k*dx2 - ClumpCenterZ;
      x2r = x2l + dx2;
      for (j=0; j<GridDimension[1]; j++) {
	x1l = gridx1l + j*dx1 - ClumpCenterY;
	x1r = x1l + dx1;
	for (i=0; i<GridDimension[0]; i++) {
	  idx = (k*GridDimension[1] + j)*GridDimension[0] + i;
	  x0l = gridx0l + i*dx0 - ClumpCenterX;
	  x0r = x0l + dx0;
	  d[0] = sqrt(x0l*x0l + x1l*x1l + x2l*x2l);
	  d[1] = sqrt(x0r*x0r + x1l*x1l + x2l*x2l);
	  d[2] = sqrt(x0l*x0l + x1r*x1r + x2l*x2l);
	  d[3] = sqrt(x0l*x0l + x1l*x1l + x2r*x2r);
	  d[4] = sqrt(x0r*x0r + x1r*x1r + x2l*x2l);
	  d[5] = sqrt(x0r*x0r + x1l*x1l + x2r*x2r);
	  d[6] = sqrt(x0l*x0l + x1r*x1r + x2r*x2r);
	  d[7] = sqrt(x0r*x0r + x1r*x1r + x2r*x2r);
	  Vin = 0.0;
	  for (l=0; l<8; l++) {
	    Vin += (d[l] <= ClumpRadius) ? 0.125 : 0.0;
	  }
	  Vout = 1.0 - Vin;
	  
	  // set remaining quantities based on approx. volume in/out of clump
	  BaryonField[RhoNum][idx] = (Vout*RhoConstOut + Vin*RhoConstIn)/DensityUnits;
	  BaryonField[TENum][idx]  = (Vout*TEConstOut + Vin*TEConstIn)/eUnits;
	  if (DualEnergyFormalism)
	    BaryonField[IENum][idx] = (Vout*IEConstOut + Vin*IEConstIn)/eUnits;
	  BaryonField[DeNum][idx]  = (Vout*DeConstOut + Vin*DeConstIn)/DensityUnits;
	  BaryonField[HINum][idx]  = (Vout*HIConstOut + Vin*HIConstIn)/DensityUnits;
	  BaryonField[HIINum][idx] = (Vout*HIIConstOut + Vin*HIIConstIn)/DensityUnits;
	  if ((NumChemicals == 3) || (MultiSpecies > 0)) {
	    BaryonField[HeINum][idx]   = (Vout*HeIConstOut + Vin*HeIConstIn)/DensityUnits;
	    BaryonField[HeIINum][idx]  = (Vout*HeIIConstOut + Vin*HeIIConstIn)/DensityUnits;
	    BaryonField[HeIIINum][idx] = (Vout*HeIIIConstOut + Vin*HeIIIConstIn)/DensityUnits;
	  }
	}
      }
    }
    
    // output debugging values to stdout
    if (debug) {
      
      printf("\n  Initializing constant fields using CGS values:\n");
      printf("  Outside the clump:\n");
      printf("        density = %g\n",RhoConstOut);
      printf("   total energy = %g\n",TEConstOut);
      if (DualEnergyFormalism)
	printf("   internal energy = %g\n",IEConstOut);
      printf("     x-velocity = %g\n",VxConst);
      printf("     y-velocity = %g\n",VyConst);
      printf("     z-velocity = %g\n",VzConst);
      printf("      radiation = %g\n",EgConst);
      printf("      electrons = %g\n",DeConstOut);
      printf("            nHI = %g\n",HIConstOut);
      printf("           nHII = %g\n",HIIConstOut);
      if ((NumChemicals == 3) || (MultiSpecies > 0)) {
	printf("           nHeI = %g\n",HeIConstOut);
	printf("          nHeII = %g\n",HeIIConstOut);
	printf("         nHeIII = %g\n",HeIIIConstOut);
      }
      
      printf("\n  Inside the clump:\n");
      printf("        density = %g\n",RhoConstIn);
      printf("   total energy = %g\n",TEConstIn);
      if (DualEnergyFormalism)
	printf("   internal energy = %g\n",IEConstIn);
      printf("     x-velocity = %g\n",VxConst);
      printf("     y-velocity = %g\n",VyConst);
      printf("     z-velocity = %g\n",VzConst);
      printf("      radiation = %g\n",EgConst);
      printf("      electrons = %g\n",DeConstIn);
      printf("            nHI = %g\n",HIConstIn);
      printf("           nHII = %g\n",HIIConstIn);
      if ((NumChemicals == 3) || (MultiSpecies > 0)) {
	printf("           nHeI = %g\n",HeIConstIn);
	printf("          nHeII = %g\n",HeIIConstIn);
	printf("         nHeIII = %g\n",HeIIIConstIn);
      }
    }

  } // end if NewData == TRUE

  return SUCCESS;

#else

  fprintf(stderr,"Error: RAD_HYDRO must be enabled for this test!\n");
  return FAIL;

#endif

}
