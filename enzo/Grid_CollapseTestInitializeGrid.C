/***********************************************************************
/
/  GRID CLASS (INITIALIZE THE GRID FOR A COLLAPSE TEST)
/
/  written by: Greg Bryan
/  date:       May, 1998
/  modified1:
/
/  PURPOSE:
/
/  RETURNS: FAIL or SUCCESS
/
************************************************************************/
 
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "CosmologyParameters.h"
 
int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, float *MassUnits, FLOAT Time);
int CosmologyComputeExpansionFactor(FLOAT time, FLOAT *a, FLOAT *dadt);
float gasdev();
 
static int CollapseTestParticleCount = 0;
 
static float CosmologySimulationInitialFractionHII   = 1.2e-5;
static float CosmologySimulationInitialFractionHeII  = 1.0e-14;
static float CosmologySimulationInitialFractionHeIII = 1.0e-17;
static float CosmologySimulationInitialFractionHM    = 2.0e-9;
static float CosmologySimulationInitialFractionH2I   = 2.0e-20;
static float CosmologySimulationInitialFractionH2II  = 3.0e-14;
 
int grid::CollapseTestInitializeGrid(int NumberOfSpheres,
			     FLOAT SphereRadius[MAX_SPHERES],
			     FLOAT SphereCoreRadius[MAX_SPHERES],
			     float SphereDensity[MAX_SPHERES],
			     float SphereTemperature[MAX_SPHERES],
			     FLOAT SpherePosition[MAX_SPHERES][MAX_DIMENSION],
			     float SphereVelocity[MAX_SPHERES][MAX_DIMENSION],
			     int   SphereType[MAX_SPHERES],
			     int   SphereUseParticles,
			     float UniformVelocity[MAX_DIMENSION],
			     int   SphereUseColour,
			     float InitialTemperature, int level)
{
  /* declarations */
 
  int dim, i, j, k, m, field, sphere, size;
  int DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum, HMNum, H2INum, H2IINum,
      DINum, DIINum, HDINum;
 
  /* create fields */
 
  NumberOfBaryonFields = 0;
  FieldType[NumberOfBaryonFields++] = Density;
  FieldType[NumberOfBaryonFields++] = TotalEnergy;
  if (DualEnergyFormalism)
    FieldType[NumberOfBaryonFields++] = InternalEnergy;
  int ivel = NumberOfBaryonFields;
  FieldType[NumberOfBaryonFields++] = Velocity1;
  if (GridRank > 1)
    FieldType[NumberOfBaryonFields++] = Velocity2;
  if (GridRank > 2)
    FieldType[NumberOfBaryonFields++] = Velocity3;
  if (MultiSpecies) {
    FieldType[DeNum    = NumberOfBaryonFields++] = ElectronDensity;
    FieldType[HINum    = NumberOfBaryonFields++] = HIDensity;
    FieldType[HIINum   = NumberOfBaryonFields++] = HIIDensity;
    FieldType[HeINum   = NumberOfBaryonFields++] = HeIDensity;
    FieldType[HeIINum  = NumberOfBaryonFields++] = HeIIDensity;
    FieldType[HeIIINum = NumberOfBaryonFields++] = HeIIIDensity;
    if (MultiSpecies > 1) {
      FieldType[HMNum    = NumberOfBaryonFields++] = HMDensity;
      FieldType[H2INum   = NumberOfBaryonFields++] = H2IDensity;
      FieldType[H2IINum  = NumberOfBaryonFields++] = H2IIDensity;
    }
    if (MultiSpecies > 2) {
      FieldType[DINum   = NumberOfBaryonFields++] = DIDensity;
      FieldType[DIINum  = NumberOfBaryonFields++] = DIIDensity;
      FieldType[HDINum  = NumberOfBaryonFields++] = HDIDensity;
    }
  }
  int ColourNum = NumberOfBaryonFields;
  if (SphereUseColour)
    FieldType[NumberOfBaryonFields++] = Metallicity; /* fake it with metals */
 
  /* Return if this doesn't concern us. */
 
  if (ProcessorNumber != MyProcessorNumber) {
    NumberOfParticles = (SphereUseParticles > 0) ? 1 : 0;
    for (dim = 0; dim < GridRank; dim++)
      NumberOfParticles *= (GridEndIndex[dim] - GridStartIndex[dim] + 1);
    return SUCCESS;
  }
 
  /* Set various units. */
 
  const double Mpc = 3.0856e24, SolarMass = 1.989e33, GravConst = 6.67e-8,
               pi = 3.14159, mh = 1.67e-24, kboltz = 1.381e-16;
  float DensityUnits=1, LengthUnits=1, TemperatureUnits=1, TimeUnits=1,
    VelocityUnits=1, MassUnits=1, CriticalDensity=1, BoxLength=1, mu=0.6;
  FLOAT a, dadt, ExpansionFactor = 1;
  GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	   &TimeUnits, &VelocityUnits, &MassUnits, Time);
  if (ComovingCoordinates) {
    CosmologyComputeExpansionFactor(Time, &a, &dadt);
    ExpansionFactor = a/(1.0+InitialRedshift);
    CriticalDensity = 2.78e11*POW(HubbleConstantNow, 2); // in Msolar/Mpc^3
    BoxLength = ComovingBoxSize*ExpansionFactor/HubbleConstantNow;  // in Mpc
  }
 
  /* Compute NFW profile-info. The parameters are SphereCoreRadius (which is
     the "knee radius") and SphereDensity (the overdensity ar the knee
     radius).  The pressure is computed by integrated the hydro-static
     equation and from this comes the temperature and the dm velocity
     dispersion. */
 
#define NFW_POINTS 500
  float NFWMass[NFW_POINTS], NFWPressure[NFW_POINTS], NFWTemp[NFW_POINTS], x1,
        NFWDensity[NFW_POINTS], NFWSigma[NFW_POINTS], m200;
  FLOAT NFWRadius[NFW_POINTS];
  double dpdr = 0, dpdr_old;
  sphere = 0;
  m200   = 0;
  NFWPressure[0] = 1.0 * kboltz * InitialTemperature / (mu * mh);
  FILE *fptr = fopen("NFWProfile.out", "w");
  for (i = 0; i < NFW_POINTS; i++) {
    NFWRadius[i] = SphereRadius[sphere]*POW(10, -3*(float(i)/NFW_POINTS));
    x1 = NFWRadius[i]/SphereCoreRadius[sphere];
    NFWDensity[i] = SphereDensity[sphere]/(x1*(1.0+x1)*(1.0+x1));
    NFWMass[i] = 4.0*pi*SphereDensity[sphere]*
                      (CriticalDensity/POW(ExpansionFactor, 3)) *
		POW(SphereCoreRadius[sphere]*BoxLength, 3) *
		(log(1.0+x1) - x1/(x1+1.0));  // in Msolar
    dpdr_old = dpdr;
    dpdr = GravConst * NFWMass[i] * SolarMass *
           NFWDensity[i] /
           POW(NFWRadius[i]*BoxLength*Mpc, 2);
    if (i > 0)
      NFWPressure[i] = NFWPressure[i-1] -
	0.5*(dpdr+dpdr_old)*(NFWRadius[i]-NFWRadius[i-1])*BoxLength*Mpc;
    NFWTemp[i] = NFWPressure[i]*mu*mh/(kboltz*NFWDensity[i]); // in K
    NFWSigma[i] = sqrt(kboltz * NFWTemp[i] / (mu * mh));  // in cm/s
    float mean_overdensity = 3.0*SphereDensity[sphere] / (x1*x1*x1) *
        (log(1.0+x1) - x1/(x1+1.0));
    fprintf(fptr, "%"ISYM" %"GOUTSYM" %"GSYM" %"GSYM" %"GSYM" %"GSYM" %"GSYM" %"GSYM"\n", i, NFWRadius[i],
	 NFWDensity[i], NFWMass[i], NFWPressure[i], NFWTemp[i], NFWSigma[i],
         mean_overdensity);
    if (mean_overdensity > 200 && m200 == 0)
      m200 = NFWMass[i];
  }
  fprintf(fptr, "#m200 = %"GSYM"\n", m200);
  fclose(fptr);
 
  /* Loop over the set-up twice, once to count the particles, the second
     time to initialize them. */
 
  int SetupLoopCount, npart = 0;
  for (SetupLoopCount = 0; SetupLoopCount < 1+min(SphereUseParticles, 1);
       SetupLoopCount++) {
 
  /* Set densities */
 
  float BaryonMeanDensity = SphereUseParticles ? 0.1 : 1.0;
  if (SphereUseParticles == 2) BaryonMeanDensity = 0.9;
  float ParticleMeanDensity = 1.0 - BaryonMeanDensity, ParticleCount = 0;
 
  /* Set particles. */
 
  if (SphereUseParticles > 0 && SetupLoopCount > 0) {
 
    /* If particles already exist (coarse particles), then delete. */
 
    if (NumberOfParticles > 0)
      this->DeleteParticles();
 
    /* Use count from previous loop to set particle number. */
 
    NumberOfParticles = npart;
    npart = 0;
 
    /* Allocate space. */
 
    this->AllocateNewParticles(NumberOfParticles);
 
    /* Particle values will be set below. */
 
  } // end: particle initialization
 
  /* Set up the baryon field. */
 
  /* compute size of fields */
 
  size = 1;
  for (dim = 0; dim < GridRank; dim++)
    size *= GridDimension[dim];
 
  /* allocate fields */
 
  if (SetupLoopCount == 0)
    for (field = 0; field < NumberOfBaryonFields; field++)
      if (BaryonField[field] == NULL)
	BaryonField[field] = new float[size];
 
  /* Loop over the mesh. */
 
  float density, dens1, Velocity[MAX_DIMENSION],
    temperature, temp1, sigma, sigma1, colour;
  FLOAT r, x, y = 0, z = 0;
  int n = 0;
 
  for (k = 0; k < GridDimension[2]; k++)
    for (j = 0; j < GridDimension[1]; j++)
      for (i = 0; i < GridDimension[0]; i++, n++) {
 
	/* Compute position */
 
	x = CellLeftEdge[0][i] + 0.5*CellWidth[0][i];
	if (GridRank > 1)
	  y = CellLeftEdge[1][j] + 0.5*CellWidth[1][j];
	if (GridRank > 2)
	  z = CellLeftEdge[2][k] + 0.5*CellWidth[2][k];
 
	/* Loop over spheres. */
 
	density = 1.0;
	temperature = temp1 = InitialTemperature;
	sigma = sigma1 = 0;
	colour = 1.0e-10;
	for (dim = 0; dim < MAX_DIMENSION; dim++)
	  Velocity[dim] = 0;
	for (sphere = 0; sphere < NumberOfSpheres; sphere++) {
 
	  /* Find distance from center. */
 
	  r = sqrt(POW(fabs(x-SpherePosition[sphere][0]), 2) +
		   POW(fabs(y-SpherePosition[sphere][1]), 2) +
		   POW(fabs(z-SpherePosition[sphere][2]), 2) );
	  r = max(r, 0.1*CellWidth[0][0]);
 
	  if (r < SphereRadius[sphere]) {
 
	    /* 1) Uniform */
 
	    if (SphereType[sphere] == 1)
	      dens1 = SphereDensity[sphere];
 
	    /* 2) r^-2 POWer law */
 
	    if (SphereType[sphere] == 2)
	      dens1 = SphereDensity[sphere]*POW(r/SphereRadius[sphere], -2);
 
	    /* 3) NFW profile (use look-up table for temperature and
	          velocity dispersion)*/
 
	    if (SphereType[sphere] == 3) {
	      x1 = r/SphereCoreRadius[sphere];
	      dens1 = SphereDensity[sphere]/(x1*(1.0+x1)*(1.0+x1));
	      for (m = 1; m < NFW_POINTS; m++)
		if (r > NFWRadius[m]) {
		  temp1 = NFWTemp[m] + (NFWTemp[m-1] - NFWTemp[m])*
		    (r - NFWRadius[m])/(NFWRadius[m-1] - NFWRadius[m]);
		  sigma1 = NFWSigma[m] + (NFWSigma[m-1] - NFWSigma[m])*
		    (r - NFWRadius[m])/(NFWRadius[m-1] - NFWRadius[m]);
		  break;
		}
	    }
 
	    /* 4) Gaussian */
 
	    if (SphereType[sphere] == 4) {
	      dens1 = SphereDensity[sphere]*
                      exp(-0.5*POW(r/SphereCoreRadius[sphere], 2));
	    }
 
	    /* 5) r^-2 POWer law with core radius */
 
	    if (SphereType[sphere] == 5) {
	      if (r < SphereCoreRadius[sphere])
		dens1 = SphereDensity[sphere]*POW(SphereCoreRadius[sphere]/
						  SphereRadius[sphere], -2);
	      else
		dens1 = SphereDensity[sphere]*POW(r/SphereRadius[sphere], -2);
	    }
 
	    /* 10) disk (ok, it's not a sphere, so shoot me) */
 
	    if (SphereType[sphere] == 10) {
 
	      FLOAT xpos, ypos, zpos, xpos1, ypos1, zpos1, zheight, drad;
	      FLOAT ScaleHeightz = SphereCoreRadius[sphere]/6.0,
	            ScaleHeightR = SphereCoreRadius[sphere];
 
	      /* Loop over dims if using Zeus (since vel's face-centered). */
 
	      for (dim = 0; dim < 1+(HydroMethod == Zeus_Hydro ? GridRank : 0);
		   dim++) {
 
		/* Compute position. */
 
		xpos = x-SpherePosition[sphere][0] -
		  (dim == 1 ? 0.5*CellWidth[0][0] : 0.0);
		ypos = y-SpherePosition[sphere][1] -
		  (dim == 2 ? 0.5*CellWidth[1][0] : 0.0);
		zpos = z-SpherePosition[sphere][2] -
		  (dim == 3 ? 0.5*CellWidth[2][0] : 0.0);
 
		/* Compute z and r_perp (SphereVelocity is angular momentum
		   and must have unit length). */
 
		zheight = SphereVelocity[sphere][0]*xpos +
		          SphereVelocity[sphere][1]*ypos +
		          SphereVelocity[sphere][2]*zpos;
		xpos1 = xpos - zheight*SphereVelocity[sphere][0];
		ypos1 = ypos - zheight*SphereVelocity[sphere][1];
		zpos1 = zpos - zheight*SphereVelocity[sphere][2];
		drad = sqrt(xpos1*xpos1 + ypos1*ypos1 + zpos1*zpos1);
 
		/* If we're above the disk, then exit. */
 
//		if (zheight > max(5.0*ScaleHeightz, 2.0*CellWidth[0][0]))
//		  continue;
 
		/* Compute density (Kruit & Searle 1982). */
 
		if (dim == 0)
		  dens1 = SphereDensity[sphere]*exp(-drad/ScaleHeightR)/
		    POW(cosh(zheight/max(ScaleHeightz, CellWidth[0][0])), 2);
 
		if (dens1 < density)
		  break;
 
		/* Compute velocity magnitude (divided by drad).
		   This assumes PointSourceGravityPosition and Sphere center
		   are the same.  This should be fixed to use the disk mass
		   as well, but that's a bit tricky. */
 
//		float vel = sqrt(PointSourceGravityConstant/drad)/drad;
 
		float accel = 0;
		if (PointSourceGravity == 1)
		  accel = PointSourceGravityConstant/
		    (POW(drad,3) + POW(PointSourceGravityCoreRadius, 3));
		if (PointSourceGravity == 2) {
		  x1 = drad/PointSourceGravityCoreRadius;
		  accel = PointSourceGravityConstant*(log(1+x1)-x1/(1+x1))/
		           POW(drad, 3);
		}
		
		float vel = sqrt(accel);
		
		/* Compute velocty: L x r_perp. */
 
		if (dim == 0 || dim == 1)
		  Velocity[0] = vel*(SphereVelocity[sphere][1]*zpos1 -
				     SphereVelocity[sphere][2]*ypos1);
		if (dim == 0 || dim == 2)
		  Velocity[1] = vel*(SphereVelocity[sphere][2]*xpos1 +
				     SphereVelocity[sphere][0]*zpos1);
		if (dim == 0 || dim == 3)
		  Velocity[2] = vel*(SphereVelocity[sphere][0]*ypos1 -
				     SphereVelocity[sphere][1]*xpos1);
 
	      } // end: loop over dims
 
	    } // end: disk
	
	    /* If the density is larger than the background (or the previous
	       sphere), then set the velocity. */
 
	    if (dens1 > density) {
	      density = dens1;
	      if (temp1 == InitialTemperature)
		temp1 = SphereTemperature[sphere];
	      temperature = temp1;
	      sigma = sigma1;
	      if (SphereType[sphere] != 10)
		for (dim = 0; dim < GridRank; dim++)
		  Velocity[dim] = SphereVelocity[sphere][dim];
	      if (sphere == 0)
		colour = dens1; /* only mark first sphere */
	    }
 
	  } // end: if (r < SphereRadius)
	} // end: loop over spheres
 
	/* Set density. */
 
	BaryonField[0][n] = density*BaryonMeanDensity;
 
	/* If doing multi-species (HI, etc.), set these. */
 
	if (MultiSpecies > 0) {
	
	  BaryonField[HIINum][n] = CosmologySimulationInitialFractionHII *
	    CoolData.HydrogenFractionByMass * BaryonField[0][n] *
	    sqrt(OmegaMatterNow)/
	    (OmegaMatterNow*BaryonMeanDensity*HubbleConstantNow);
	    //	    (CosmologySimulationOmegaBaryonNow*HubbleConstantNow);
 
	  BaryonField[HeIINum][n] = CosmologySimulationInitialFractionHeII*
	    BaryonField[0][n] * 4.0 * (1.0-CoolData.HydrogenFractionByMass);
	  BaryonField[HeIIINum][n] = CosmologySimulationInitialFractionHeIII*
	    BaryonField[0][n] * 4.0 * (1.0-CoolData.HydrogenFractionByMass);
	  BaryonField[HeINum][n] =
	    (1.0 - CoolData.HydrogenFractionByMass)*BaryonField[0][n] -
	    BaryonField[HeIINum][n] - BaryonField[HeIIINum][n];
 
	  if (MultiSpecies > 1) {
	    BaryonField[HMNum][n] = CosmologySimulationInitialFractionHM*
	      BaryonField[HIINum][n]* POW(temperature,float(0.88));
	    BaryonField[H2IINum][n] = CosmologySimulationInitialFractionH2II*
	      2.0*BaryonField[HIINum][n]* POW(temperature,float(1.8));
	    BaryonField[H2INum][n] = CosmologySimulationInitialFractionH2I*
	      BaryonField[0][n]*CoolData.HydrogenFractionByMass*POW(301.0,5.1)*
	      POW(OmegaMatterNow, float(1.5))/
	      (OmegaMatterNow*BaryonMeanDensity)/
	    //	      CosmologySimulationOmegaBaryonNow/
	      HubbleConstantNow*2.0;
	  }
 
	  BaryonField[HINum][n] =
	    CoolData.HydrogenFractionByMass*BaryonField[0][n]
	    - BaryonField[HIINum][n];
	  if (MultiSpecies > 1)
	    BaryonField[HINum][n] -= BaryonField[HMNum][n]
	      + BaryonField[H2IINum][n]
	      + BaryonField[H2INum][n];
 
	  BaryonField[DeNum][n] = BaryonField[HIINum][n] +
	    0.25*BaryonField[HeIINum][n] + 0.5*BaryonField[HeIIINum][n];
	  if (MultiSpecies > 1)
	    BaryonField[DeNum][n] += 0.5*BaryonField[H2IINum][n] -
	      BaryonField[HMNum][n];
 
	  /* Set Deuterium species (assumed to be negligible). */
 
	  if (MultiSpecies > 2) {
	    BaryonField[DINum][n] = CoolData.DeuteriumToHydrogenRatio*
	                              BaryonField[HINum][n];
	    BaryonField[DIINum][n] = CoolData.DeuteriumToHydrogenRatio*
	                             BaryonField[HIINum][n];
	    BaryonField[HDINum][n] = CoolData.DeuteriumToHydrogenRatio*
	                             BaryonField[H2INum][n];
	  }
	}
 
	/* If there is a colour field, set it. */
 
	if (SphereUseColour)
	  BaryonField[ColourNum][n] = colour;
 
	/* Set Velocities. */
 
	for (dim = 0; dim < GridRank; dim++)
	  BaryonField[ivel+dim][n] = Velocity[dim] + UniformVelocity[dim];
 
	/* Set energy (thermal and then total if necessary). */
 
	BaryonField[1][n] = temperature/TemperatureUnits/
                            ((Gamma-1.0)*mu);
 
	if (DualEnergyFormalism)
	  BaryonField[2][n] = BaryonField[1][n];
 
	if (HydroMethod != Zeus_Hydro)
	  for (dim = 0; dim < GridRank; dim++)
	    BaryonField[1][n] += 0.5*POW(BaryonField[ivel+dim][n], 2);
 
	/* Set particles if being used (generate a number of particle
	   proportional to density). */
 
	if (SphereUseParticles)
	  if (i >= GridStartIndex[0] && i <= GridEndIndex[0] &&
	      j >= GridStartIndex[1] && j <= GridEndIndex[1] &&
	      k >= GridStartIndex[2] && k <= GridEndIndex[2]  ) {
	    ParticleCount += density/POW(float(RefineBy), GridRank*level);
	    while (ParticleCount > 1) {
	      if (SetupLoopCount > 0) {
		ParticleMass[npart] = ParticleMeanDensity*
		                      POW(float(RefineBy), GridRank*level);
		ParticleNumber[npart] = CollapseTestParticleCount++;
		ParticleType[npart] = PARTICLE_TYPE_DARK_MATTER;
 
		/* Set random position within cell. */
 
		ParticlePosition[0][npart] = x +
		        CellWidth[0][0]*(FLOAT(rand())/FLOAT(RAND_MAX) - 0.5);
		ParticlePosition[1][npart] = y +
		        CellWidth[1][0]*(FLOAT(rand())/FLOAT(RAND_MAX) - 0.5);
		ParticlePosition[2][npart] = z +
		        CellWidth[2][0]*(FLOAT(rand())/FLOAT(RAND_MAX) - 0.5);
 
		/* Set bulk velocity. */
 
		for (dim = 0; dim < GridRank; dim++)
		  ParticleVelocity[dim][npart] =
		    Velocity[dim]+UniformVelocity[dim];
 
		/* Add random velocity; */
 
		if (sigma != 0)
		  for (dim = 0; dim < GridRank; dim++)
		    ParticleVelocity[dim][npart] +=
		      gasdev()*sigma/VelocityUnits;
 
	      }
	      npart++;
	      ParticleCount -= 1.0;
	    }
	  } // end: if statement
 
      } // end loop over grid
 
  } // end loop SetupLoopCount
 
  if (SphereUseParticles && debug)
    printf("CollapseTestInitialize: NumberOfParticles = %"ISYM"\n",
	   NumberOfParticles);
 
  return SUCCESS;
}
 
 
/* Routine to return a Gaussian random deviate with zero mean and unite
   variance (adapted from Numerical Recipes). */
 
static int gasdev_iset = 0;
static float gasdev_gset;
 
float gasdev()
{
  float v1, v2, r = 0, fac, gasdev_ret;
  if (gasdev_iset == 0) {
    while (r >= 1 || r == 0) {
      v1 = 2.0*float(rand())/(float(RAND_MAX)) - 1.0;
      v2 = 2.0*float(rand())/(float(RAND_MAX)) - 1.0;
      r = v1*v1 + v2*v2;
    }
    fac = sqrt(-2.0*log(r)/r);
    gasdev_gset = v1*fac;
    gasdev_ret  = v2*fac;
    gasdev_iset = 1;
  } else {
    gasdev_ret  = gasdev_gset;
    gasdev_iset = 0;
  }
  return gasdev_ret;
}
