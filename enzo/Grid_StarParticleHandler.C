/***********************************************************************
/
/  GRID CLASS (HANDLE THE CREATION AND FEEDBACK OF STAR PARTICLES)
/
/  written by: Greg Bryan
/  date:       March, 1997
/  modified1:
/
/  PURPOSE:
/
/  RETURNS:
/    SUCCESS or FAIL
/
************************************************************************/
 
#include <stdio.h>
#include <math.h>
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "fortran.def"
#include "CosmologyParameters.h"
#include "StarParticleData.h"

#define  PROTONMASS  1.6726e-24

/* function prototypes */
 
int CosmologyComputeExpansionFactor(FLOAT time, FLOAT *a, FLOAT *dadt);
int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, float *MassUnits, FLOAT Time);
int FindField(int field, int farray[], int numfields);
 
#define NO_STAR1
 
#ifdef STAR1
extern "C" void FORTRAN_NAME(star_maker1)(int *nx, int *ny, int *nz,
             float *d, float *dm, float *temp, float *u, float *v, float *w,
                float *cooltime,
             float *dt, float *r, float *dx, FLOAT *t, float *z, int *procnum,
             float *d1, float *x1, float *v1, float *t1,
             int *nmax, FLOAT *xstart, FLOAT *ystart, FLOAT *zstart,
     		 int *ibuff, hydro_method *imethod,
             float *odthresh, float *massff, float *smthrest, int *level,
                 int *np,
             FLOAT *xp, FLOAT *yp, FLOAT *zp, float *up, float *vp, float *wp,
             float *mp, float *tdp, float *tcp);
#endif /* STAR1 */
 
extern "C" void FORTRAN_NAME(star_maker2)(int *nx, int *ny, int *nz,
             float *d, float *dm, float *temp, float *u, float *v, float *w,
                float *cooltime,
             float *dt, float *r, float *metal, float *dx, FLOAT *t, float *z,
             int *procnum,
             float *d1, float *x1, float *v1, float *t1,
             int *nmax, FLOAT *xstart, FLOAT *ystart, FLOAT *zstart,
     		 int *ibuff,
             int *imetal, hydro_method *imethod, float *mintdyn,
             float *odthresh, float *massff, float *smthrest, int *level,
                 int *np,
             FLOAT *xp, FLOAT *yp, FLOAT *zp, float *up, float *vp, float *wp,
             float *mp, float *tdp, float *tcp, float *metalf);
 
extern "C" void FORTRAN_NAME(star_maker3)(int *nx, int *ny, int *nz,
             float *d, float *dm, float *temp, float *u, float *v, float *w,
                float *cooltime,
             float *dt, float *r, float *metal, float *zfield1, float *zfield2,
             float *dx, FLOAT *t, float *z,
             int *procnum,
             float *d1, float *x1, float *v1, float *t1,
             int *nmax, FLOAT *xstart, FLOAT *ystart, FLOAT *zstart,
     		 int *ibuff,
             int *imetal, hydro_method *imethod, float *mintdyn,
             float *odthresh, float *massff, float *smthrest, int *level,
                 int *np,
             FLOAT *xp, FLOAT *yp, FLOAT *zp, float *up, float *vp, float *wp,
             float *mp, float *tdp, float *tcp, float *metalf);
 
extern "C" void FORTRAN_NAME(star_maker4)(int *nx, int *ny, int *nz,
             float *d, float *dm, float *temp, float *u, float *v, float *w,
                float *cooltime,
             float *dt, float *r, float *metal, float *dx, FLOAT *t, float *z, 
             int *procnum,
             float *d1, float *x1, float *v1, float *t1,
             int *nmax, FLOAT *xstart, FLOAT *ystart, FLOAT *zstart, 
     		 int *ibuff, 
             int *imetal, hydro_method *imethod, float *mintdyn,
             float *odthresh, float *massff, float *smthrest, int *level,
                 int *np,
             FLOAT *xp, FLOAT *yp, FLOAT *zp, float *up, float *vp, float *wp,
             float *mp, float *tdp, float *tcp, float *metalf);

extern "C" void FORTRAN_NAME(star_maker5)(int *nx, int *ny, int *nz,
             float *d, float *dm, float *temp, float *coolrate, float *u, float *v, float *w,
                float *cooltime,
             float *dt, float *r, float *metal, float *dx, FLOAT *t, float *z, 
             int *procnum,
             float *d1, float *x1, float *v1, float *t1,
             int *nmax, FLOAT *xstart, FLOAT *ystart, FLOAT *zstart, 
     		 int *ibuff, 
             int *imetal, hydro_method *imethod, float *mintdyn,
             float *odthresh, float *massff, float *smthrest, int *level,
                 int *np,
             FLOAT *xp, FLOAT *yp, FLOAT *zp, float *up, float *vp, float *wp,
	     float *mp, float *tdp, float *tcp, float *metalf, int ran1_init);

/* adding distributed starmaker and feedback */
extern "C" void FORTRAN_NAME(star_maker10)(int *nx, int *ny, int *nz,
             float *d, float *dm, float *temp, float *u, float *v, float *w,
                float *cooltime,
             float *dt, float *r, float *metal, float *zfield1, float *zfield2,
             float *dx, FLOAT *t, float *z,
             int *procnum,
             float *d1, float *x1, float *v1, float *t1,
             int *nmax, FLOAT *xstart, FLOAT *ystart, FLOAT *zstart,
     		 int *ibuff,
             int *imetal, hydro_method *imethod, float *mintdyn,
             float *odthresh, float *massff, float *smthrest, int *level,
		 int *np, 
             FLOAT *xp, FLOAT *yp, FLOAT *zp, float *up, float *vp, float *wp,
             float *mp, float *tdp, float *tcp, float *metalf);

#ifdef STAR1
extern "C" void FORTRAN_NAME(star_feedback1)(int *nx, int *ny, int *nz,
             float *d, float *dm, float *temp, float *u, float *v,
		       float *w, float *dt, float *r, float *dx,
                       FLOAT *t, float *z,
             float *d1, float *x1, float *v1, float *t1,
             int *nmax, FLOAT *xstart, FLOAT *ystart, FLOAT *zstart,
     				 int *ibuff,
             FLOAT *xp, FLOAT *yp, FLOAT *zp, float *up, float *vp, float *wp,
             float *mp, float *tdp, float *tcp, float *te, float *ge,
		       int *idual);
#endif /* STAR1 */
 
extern "C" void FORTRAN_NAME(star_feedback2)(int *nx, int *ny, int *nz,
             float *d, float *dm, float *te, float *ge, float *u, float *v,
		       float *w, float *metal,
             int *idual, int *imetal, hydro_method *imethod, float *dt,
		       float *r, float *dx, FLOAT *t, float *z,
             float *d1, float *x1, float *v1, float *t1,
                       float *sn_param, float *m_eject, float *yield,
             int *nmax, FLOAT *xstart, FLOAT *ystart, FLOAT *zstart,
		       int *ibuff,
             FLOAT *xp, FLOAT *yp, FLOAT *zp, float *up, float *vp, float *wp,
             float *mp, float *tdp, float *tcp, float *metalf,
			float *justburn);
 
extern "C" void FORTRAN_NAME(star_feedback3)(int *nx, int *ny, int *nz,
             float *d, float *dm, float *te, float *ge, float *u, float *v,
		       float *w, float *metal, float *zfield1, float *zfield2,
	     int *idual, int *imetal, int *imulti_metals, hydro_method *imethod, 
		       float *dt, float *r, float *dx, FLOAT *t, float *z,
             float *d1, float *x1, float *v1, float *t1,
                       float *sn_param, float *m_eject, float *yield,
             int *nmax, FLOAT *xstart, FLOAT *ystart, FLOAT *zstart,
		       int *ibuff,
             FLOAT *xp, FLOAT *yp, FLOAT *zp, float *up, float *vp, float *wp,
             float *mp, float *tdp, float *tcp, float *metalf,
			float *justburn);
 
extern "C" void FORTRAN_NAME(star_feedback4)(int *nx, int *ny, int *nz,
             float *d, float *dm, float *te, float *ge, float *u, float *v, 
		       float *w, float *metal,
             int *idual, int *imetal, hydro_method *imethod, float *dt, 
		       float *r, float *dx, FLOAT *t, float *z, 
             float *d1, float *x1, float *v1, float *t1,
                       float *sn_param, float *m_eject, float *yield,
             int *nmax, FLOAT *xstart, FLOAT *ystart, FLOAT *zstart, 
		       int *ibuff,
             FLOAT *xp, FLOAT *yp, FLOAT *zp, float *up, float *vp, float *wp,
             float *mp, float *tdp, float *tcp, float *metalf, 
			float *justburn);

extern "C" void FORTRAN_NAME(star_feedback5)(int *nx, int *ny, int *nz,
             float *d, float *dm, float *te, float *ge, float *u, float *v, 
		       float *w, float *metal,
             int *idual, int *imetal, hydro_method *imethod, float *dt, 
		       float *r, float *dx, FLOAT *t, float *z, 
             float *d1, float *x1, float *v1, float *t1,
                       float *sn_param, float *m_eject, float *yield,
             int *nmax, FLOAT *xstart, FLOAT *ystart, FLOAT *zstart, 
		       int *ibuff,
             FLOAT *xp, FLOAT *yp, FLOAT *zp, float *up, float *vp, float *wp,
	     float *mp, float *tdp, float *tcp, float *metalf,
			float *justburn);
 
extern "C" void FORTRAN_NAME(copy3d)(float *source, float *dest,
                                   int *sdim1, int *sdim2, int *sdim3,
                                   int *ddim1, int *ddim2, int *ddim3,
                                   int *sstart1, int *sstart2, int *sstart3,
                                   int *dstart1, int *dstart2, int *dststart3);

/* adding distributed starmaker and feedback */
extern "C" void FORTRAN_NAME(star_feedback10)(int *nx, int *ny, int *nz,
             float *d, float *dm, float *te, float *ge, float *u, float *v,
		       float *w, float *metal, float *zfield1, float *zfield2,
	     int *idual, int *imetal, int *imulti_metals, hydro_method *imethod, 
		       float *dt, float *r, float *dx, FLOAT *t, float *z,
             float *d1, float *x1, float *v1, float *t1,
                       float *sn_param, float *m_eject, float *yield,
	     int *distrad, int *diststep, int *distcells,
             int *nmax, FLOAT *xstart, FLOAT *ystart, FLOAT *zstart,
		       int *ibuff,
             FLOAT *xp, FLOAT *yp, FLOAT *zp, float *up, float *vp, float *wp,
             float *mp, float *tdp, float *tcp, float *metalf, int *type,
			float *justburn);
 
#ifdef EMISSIVITY
  int CalcEmiss(int *nx, int *ny, int *nz,
             float *d, float *dm, float *te, float *ge, float *u, float *v,
		       float *w, float *metal,
             int *idual, int *imetal, hydro_method *imethod, float *dt,
		       float *r, float *dx, FLOAT *t, float *z,
             float *d1, float *x1, float *v1, float *t1,
                       float *sn_param, float *m_eject, float *yield,
             int *nmax, FLOAT *xstart, FLOAT *ystart, FLOAT *zstart,
		       int *ibuff,
             FLOAT *xp, FLOAT *yp, FLOAT *zp, float *up, float *vp, float *wp,
             float *mp, float *tdp, float *tcp, float *metalf,
	      float *justburn, float *EmissivityArray, float dtLevelAbove);
#endif
 
 
#ifdef EMISSIVITY
int grid::StarParticleHandler(int level, float dtLevelAbove)
#else
int grid::StarParticleHandler(int level)
#endif
{
 
  if (MyProcessorNumber != ProcessorNumber)
    return SUCCESS;
 
  if (NumberOfBaryonFields == 0)
    return SUCCESS;
 
  /* initialize */
 
  int dim, i, j, k, index, size, field, GhostZones = DEFAULT_GHOST_ZONES;
  int DensNum, GENum, TENum, Vel1Num, Vel2Num, Vel3Num;
 
  /* Compute size (in floats) of the current grid. */
 
  size = 1;
  for (dim = 0; dim < GridRank; dim++)
    size *= GridDimension[dim];
 
  /* Find fields: density, total energy, velocity1-3. */
 
  this->DebugCheck("StarParticleHandler");
  if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num,
				       Vel3Num, TENum) == FAIL) {
    fprintf(stderr, "Error in IdentifyPhysicalQuantities.\n");
    return FAIL;
  }
 
  /* Find metallicity field and set flag. */
 
  int MetallicityField = FALSE, MetalNum;
  if ((MetalNum = FindField(Metallicity, FieldType, NumberOfBaryonFields))
      != -1)
    MetallicityField = TRUE;
  else
    MetalNum = 0;
 
  /* Compute the redshift. */
 
  float zred;
  FLOAT a = 1, dadt;
  if (ComovingCoordinates)
    CosmologyComputeExpansionFactor(Time, &a, &dadt);
  zred = 1.0*(1.0+InitialRedshift)/a - 1.0;
 
  /* Compute the temperature field. */
 
  float *temperature = new float[size];
  this->ComputeTemperatureField(temperature);
 
  /* Get the dark matter field in a usable size for star_maker
     (if level > MaximumGravityRefinementLevel then the dark matter
      field is not valid, so just make it zero - by this time, the
      evolution will be dominated by baryonic matter anyway). */
 
  float *dmfield = new float[size];
  int StartIndex[MAX_DIMENSION], Zero[] = {0,0,0};
  if (level <= MaximumGravityRefinementLevel &&
      GravitatingMassFieldParticles != NULL) {
    for (dim = 0; dim < MAX_DIMENSION; dim++)
      StartIndex[dim] =
      nint((CellLeftEdge[dim][0] - GravitatingMassFieldParticlesLeftEdge[dim])/
	   GravitatingMassFieldParticlesCellSize);
    FORTRAN_NAME(copy3d)(GravitatingMassFieldParticles, dmfield,
			 GravitatingMassFieldParticlesDimension,
			 GravitatingMassFieldParticlesDimension+1,
			 GravitatingMassFieldParticlesDimension+2,
			 GridDimension, GridDimension+1, GridDimension+2,
			 Zero, Zero+1, Zero+2,
			 StartIndex, StartIndex+1, StartIndex+2);
  } else
    for (i = 0; i < size; i++)
      dmfield[i] = 0;
 
  /* Convert the species densities into fractional densities (i.e. divide
     by total baryonic density).  At the end we will multiply by the new
     density so that species fractions are maintained. */
 
  for (field = 0; field < NumberOfBaryonFields; field++) {

    if (FieldType[field] >= ElectronDensity &&
#ifdef EMISSIVITY
	FieldType[field] <  GravPotential
#else
	FieldType[field] <  FieldUndefined
#endif
        ) {

      for (k = GridStartIndex[2]; k <= GridEndIndex[2]; k++) {
	for (j = GridStartIndex[1]; j <= GridEndIndex[1]; j++) {
	  index = (k*GridDimension[1] + j)*GridDimension[0] +
	    GridStartIndex[0];
	  for (i = GridStartIndex[0]; i <= GridEndIndex[0]; i++, index++)
	    BaryonField[field][index] /= BaryonField[DensNum][index];
	}
      }

    } // if Density

  }
 
  /* Set the units. */
 
  float DensityUnits = 1, LengthUnits = 1, TemperatureUnits = 1,
    TimeUnits = 1, VelocityUnits = 1, MassUnits = 1;
  if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	       &TimeUnits, &VelocityUnits, &MassUnits, Time) == FAIL) {
    fprintf(stderr, "Error in GetUnits.\n");
    return FAIL;    
  }
 
  float CellWidthTemp = float(CellWidth[0][0]);
 
  /* ------------------------------------------------------------------- */
  /* 1) StarParticle creation. */
 
  //  if (StarParticleCreation > 0 && level == MaximumRefinementLevel) {
  if (StarParticleCreation > 0) {
 
    /* Generate a fake grid to keep the particles in. */
 
    grid *tg = new grid;
    tg->GridRank = GridRank;
    tg->ProcessorNumber = ProcessorNumber;
 
    /* Allocate space for new particles. */
 
    int MaximumNumberOfNewParticles = int(0.25*float(size)) + 5;
    tg->AllocateNewParticles(MaximumNumberOfNewParticles);
 
    /* Compute the cooling time. */
 
    float *cooling_time = new float[size];
    this->ComputeCoolingTime(cooling_time);
 
    /* Call FORTRAN routine to do the actual work. */
 
    int NumberOfNewParticles = 0;
 
    if (debug) {
       fprintf(stderr, "StarParticle: Before Call: New StarParticles = "
	       "%"ISYM"\n", NumberOfNewParticles);
    }

#ifdef STAR1
    if (StarParticleCreation == 1) {
      FORTRAN_NAME(star_maker1)(
       GridDimension, GridDimension+1, GridDimension+2,
          BaryonField[DensNum], dmfield,
          temperature, BaryonField[Vel1Num],
          BaryonField[Vel2Num], BaryonField[Vel3Num], cooling_time,
       &dtFixed, BaryonField[NumberOfBaryonFields], &CellWidthTemp,
          &Time, &zred, &MyProcessorNumber,
       &DensityUnits, &LengthUnits, &VelocityUnits, &TimeUnits,
       &MaximumNumberOfNewParticles, CellLeftEdge[0], CellLeftEdge[1],
          CellLeftEdge[2], &GhostZones, &HydroMethod,
       &StarMakerOverDensityThreshold, &StarMakerMassEfficiency,
          &StarMakerMinimumMass, &level, &NumberOfNewParticles,
       tg->ParticlePosition[0], tg->ParticlePosition[1],
          tg->ParticlePosition[2],
       tg->ParticleVelocity[0], tg->ParticleVelocity[1],
          tg->ParticleVelocity[2],
       tg->ParticleMass, tg->ParticleAttribute[1], tg->ParticleAttribute[0]);
    }
#endif /* STAR1 */
 
    if (StarParticleCreation == 2) {

      //---- MODIFIED SF ALGORITHM ("STANDARD VERSION")

      FORTRAN_NAME(star_maker2)(
       GridDimension, GridDimension+1, GridDimension+2,
       BaryonField[DensNum], dmfield, temperature, BaryonField[Vel1Num],
          BaryonField[Vel2Num], BaryonField[Vel3Num], cooling_time,
       &dtFixed, BaryonField[NumberOfBaryonFields], BaryonField[MetalNum],
          &CellWidthTemp, &Time, &zred, &MyProcessorNumber,
       &DensityUnits, &LengthUnits, &VelocityUnits, &TimeUnits,
       &MaximumNumberOfNewParticles, CellLeftEdge[0], CellLeftEdge[1],
          CellLeftEdge[2], &GhostZones,
       &MetallicityField, &HydroMethod, &StarMakerMinimumDynamicalTime,
       &StarMakerOverDensityThreshold, &StarMakerMassEfficiency,
          &StarMakerMinimumMass, &level, &NumberOfNewParticles,
       tg->ParticlePosition[0], tg->ParticlePosition[1],
          tg->ParticlePosition[2],
       tg->ParticleVelocity[0], tg->ParticleVelocity[1],
          tg->ParticleVelocity[2],
       tg->ParticleMass, tg->ParticleAttribute[1], tg->ParticleAttribute[0],
          tg->ParticleAttribute[2]);
    } 
    if (StarParticleCreation == 3) {

      //---- UNIGRID ALGORITHM (NO JEANS MASS)

      FORTRAN_NAME(star_maker3)(
       GridDimension, GridDimension+1, GridDimension+2,
       BaryonField[DensNum], dmfield, temperature, BaryonField[Vel1Num],
          BaryonField[Vel2Num], BaryonField[Vel3Num], cooling_time,
       &dtFixed, BaryonField[NumberOfBaryonFields], BaryonField[MetalNum],
       BaryonField[MetalNum+1], BaryonField[MetalNum+2],
          &CellWidthTemp, &Time, &zred, &MyProcessorNumber,
       &DensityUnits, &LengthUnits, &VelocityUnits, &TimeUnits,
       &MaximumNumberOfNewParticles, CellLeftEdge[0], CellLeftEdge[1],
          CellLeftEdge[2], &GhostZones,
       &MetallicityField, &HydroMethod, &StarMakerMinimumDynamicalTime,
       &StarMakerOverDensityThreshold, &StarMakerMassEfficiency,
          &StarMakerMinimumMass, &level, &NumberOfNewParticles,
       tg->ParticlePosition[0], tg->ParticlePosition[1],
          tg->ParticlePosition[2],
       tg->ParticleVelocity[0], tg->ParticleVelocity[1],
          tg->ParticleVelocity[2],
       tg->ParticleMass, tg->ParticleAttribute[1], tg->ParticleAttribute[0],
          tg->ParticleAttribute[2]);
    }

    if (StarParticleCreation == 4) {

      //---- KRAVTSOV SF ALGORITHM

      FORTRAN_NAME(star_maker4)(
       GridDimension, GridDimension+1, GridDimension+2,
       BaryonField[DensNum], dmfield, temperature, BaryonField[Vel1Num],
          BaryonField[Vel2Num], BaryonField[Vel3Num], cooling_time,
       &dtFixed, BaryonField[NumberOfBaryonFields], BaryonField[MetalNum], 
          &CellWidthTemp, &Time, &zred, &MyProcessorNumber,
       &DensityUnits, &LengthUnits, &VelocityUnits, &TimeUnits,
       &MaximumNumberOfNewParticles, CellLeftEdge[0], CellLeftEdge[1],
          CellLeftEdge[2], &GhostZones, 
       &MetallicityField, &HydroMethod, &StarMakerMinimumDynamicalTime, 
       &StarMakerOverDensityThreshold, &StarMakerMassEfficiency,
          &StarMakerMinimumMass, &level, &NumberOfNewParticles, 
       tg->ParticlePosition[0], tg->ParticlePosition[1], 
          tg->ParticlePosition[2], 
       tg->ParticleVelocity[0], tg->ParticleVelocity[1], 
          tg->ParticleVelocity[2], 
       tg->ParticleMass, tg->ParticleAttribute[1], tg->ParticleAttribute[0],
          tg->ParticleAttribute[2]);
    }


    if (StarParticleCreation == 5) {

      //---- SPRINGEL & HERNQUIST SF ALGORITHM

      float *coolingrate = new float[size];

      for(int coolindex=0; coolindex<size; coolindex++){

	float cgsdensity = BaryonField[DensNum][coolindex]*DensityUnits;
	float *electronguessptr;
	float electronguess = 0.01;
	electronguessptr = &electronguess;
	coolingrate[coolindex] = GadgetCoolingRate
	  (log10(temperature[coolindex]), cgsdensity, electronguessptr, zred);
      }						  
      
      FORTRAN_NAME(star_maker5) (
       GridDimension, GridDimension+1, GridDimension+2,
       BaryonField[DensNum], dmfield, temperature, coolingrate, BaryonField[Vel1Num],
          BaryonField[Vel2Num], BaryonField[Vel3Num], cooling_time,
       &dtFixed, BaryonField[NumberOfBaryonFields], BaryonField[MetalNum], 
          &CellWidthTemp, &Time, &zred, &MyProcessorNumber,
       &DensityUnits, &LengthUnits, &VelocityUnits, &TimeUnits,
       &MaximumNumberOfNewParticles, CellLeftEdge[0], CellLeftEdge[1],
          CellLeftEdge[2], &GhostZones, 
       &MetallicityField, &HydroMethod, &StarMakerMinimumDynamicalTime, 
       &StarMakerOverDensityThreshold, &StarMakerMassEfficiency,
          &StarMakerMinimumMass, &level, &NumberOfNewParticles, 
       tg->ParticlePosition[0], tg->ParticlePosition[1], 
          tg->ParticlePosition[2], 
       tg->ParticleVelocity[0], tg->ParticleVelocity[1], 
          tg->ParticleVelocity[2], 
       tg->ParticleMass, tg->ParticleAttribute[1], tg->ParticleAttribute[0],
       tg->ParticleAttribute[2], ran1_init);

      delete [] coolingrate;
      
    }

/* adding distributed starmaker and feedback */
    if (StarParticleCreation == 10) {

      //---- UNIGRID ALGORITHM (NO JEANS MASS)

      FORTRAN_NAME(star_maker10)(
       GridDimension, GridDimension+1, GridDimension+2,
       BaryonField[DensNum], dmfield, temperature, BaryonField[Vel1Num],
          BaryonField[Vel2Num], BaryonField[Vel3Num], cooling_time,
       &dtFixed, BaryonField[NumberOfBaryonFields], BaryonField[MetalNum],
       BaryonField[MetalNum+1], BaryonField[MetalNum+2],
          &CellWidthTemp, &Time, &zred, &MyProcessorNumber,
       &DensityUnits, &LengthUnits, &VelocityUnits, &TimeUnits,
       &MaximumNumberOfNewParticles, CellLeftEdge[0], CellLeftEdge[1],
          CellLeftEdge[2], &GhostZones,
       &MetallicityField, &HydroMethod, &StarMakerMinimumDynamicalTime,
       &StarMakerOverDensityThreshold, &StarMakerMassEfficiency,
       &StarMakerMinimumMass, &level, &NumberOfNewParticles, 
       tg->ParticlePosition[0], tg->ParticlePosition[1],
          tg->ParticlePosition[2],
       tg->ParticleVelocity[0], tg->ParticleVelocity[1],
          tg->ParticleVelocity[2],
       tg->ParticleMass, tg->ParticleAttribute[1], tg->ParticleAttribute[0],
          tg->ParticleAttribute[2]);

    }

    if (debug)
      fprintf(stderr,"StarParticle: After Formation: New StarParticles = %"ISYM"\n", NumberOfNewParticles);
  
    /* If not set in the above routine, then set the metal fraction to zero. */
 
    if (MetallicityField == FALSE || StarParticleCreation < 2)
      for (i = 0; i < NumberOfNewParticles; i++)
	tg->ParticleAttribute[2][i] = 0.0;
 
    delete [] cooling_time;
 
    /* Move any new particles into their new homes. */
 
    if (NumberOfNewParticles > 0) {
 
      if (debug)
	printf("StarParticle: New StarParticles = %"ISYM"\n", NumberOfNewParticles);
 
      /* Set the particle numbers. */
 
      for (i = 0; i < NumberOfNewParticles; i++) {
	tg->ParticleNumber[i] = INT_UNDEFINED;
	tg->ParticleType[i] = PARTICLE_TYPE_STAR;
      }

      /* Move Particles into this grid (set cell size) using the fake grid. */
 
      tg->NumberOfParticles = NumberOfNewParticles;
      for (dim = 0; dim < GridRank; dim++) {
	tg->CellWidth[dim] = new FLOAT[1];
	tg->CellWidth[dim][0] = CellWidth[dim][0];
      }
      this->MoveAllParticles(1, &tg);
 
    } // end: if (NumberOfNewParticles > 0)

#ifdef EMISSIVITY
    if (StarMakerEmissivityField > 0) {

      int EtaNum = FindField(EmissivityField0, FieldType, NumberOfBaryonFields);

      CalcEmiss(GridDimension, GridDimension+1, GridDimension+2,
          BaryonField[DensNum], dmfield,
          BaryonField[TENum], BaryonField[GENum], BaryonField[Vel1Num],
          BaryonField[Vel2Num], BaryonField[Vel3Num], BaryonField[MetalNum],
       &DualEnergyFormalism, &MetallicityField, &HydroMethod,
       &dtFixed, BaryonField[NumberOfBaryonFields], &CellWidthTemp,
          &Time, &zred,
       &DensityUnits, &LengthUnits, &VelocityUnits, &TimeUnits,
          &StarEnergyToThermalFeedback, &StarMassEjectionFraction,
          &StarMetalYield,
       &NumberOfParticles,
          CellLeftEdge[0], CellLeftEdge[1], CellLeftEdge[2], &GhostZones,
       ParticlePosition[0], ParticlePosition[1],
          ParticlePosition[2],
       ParticleVelocity[0], ParticleVelocity[1],
          ParticleVelocity[2],
       ParticleMass, ParticleAttribute[1], ParticleAttribute[0],
	  ParticleAttribute[2], &RadiationData.IntegratedStarFormation, 
       BaryonField[EtaNum], dtLevelAbove);

    }
#endif

    /* Clean up. */
 
    delete tg; // temporary grid
    this->DebugCheck("StarParticleHandler(after)");
    //    if (debug) printf("StarParticle: end\n");
 
  }
 

  /* ------------------------------------------------------------------- */
  /* 2) StarParticle feedback. */
 
#ifdef STAR1
  if (StarParticleFeedback == 1) {

    //---- THIS IS THE ORIGINAL ENZO STAR FORMATION ALG.

      FORTRAN_NAME(star_feedback1)(
       GridDimension, GridDimension+1, GridDimension+2,
          BaryonField[DensNum], dmfield, temperature,
          BaryonField[Vel1Num],
          BaryonField[Vel2Num], BaryonField[Vel3Num],
       &dtFixed, BaryonField[NumberOfBaryonFields], &CellWidthTemp,
          &Time, &zred,
       &DensityUnits, &LengthUnits, &VelocityUnits, &TimeUnits,
       &NumberOfParticles,
          CellLeftEdge[0], CellLeftEdge[1], CellLeftEdge[2], &GhostZones,
       ParticlePosition[0], ParticlePosition[1],
          ParticlePosition[2],
       ParticleVelocity[0], ParticleVelocity[1],
          ParticleVelocity[2],
       ParticleMass, ParticleAttribute[1], ParticleAttribute[0],
        BaryonField[TENum], BaryonField[GENum], &DualEnergyFormalism);
 
  } // end: if (StarParticleFeedback == 1)
#endif /* STAR1 */
 
  if (StarParticleFeedback == 2) {

    //---- THIS IS THE MODIFIED STAR FORMATION ALGORITHM
 
      FORTRAN_NAME(star_feedback2)(
       GridDimension, GridDimension+1, GridDimension+2,
          BaryonField[DensNum], dmfield,
          BaryonField[TENum], BaryonField[GENum], BaryonField[Vel1Num],
          BaryonField[Vel2Num], BaryonField[Vel3Num], BaryonField[MetalNum],
       &DualEnergyFormalism, &MetallicityField, &HydroMethod,
       &dtFixed, BaryonField[NumberOfBaryonFields], &CellWidthTemp,
          &Time, &zred,
       &DensityUnits, &LengthUnits, &VelocityUnits, &TimeUnits,
          &StarEnergyToThermalFeedback, &StarMassEjectionFraction,
          &StarMetalYield,
       &NumberOfParticles,
          CellLeftEdge[0], CellLeftEdge[1], CellLeftEdge[2], &GhostZones,
       ParticlePosition[0], ParticlePosition[1],
          ParticlePosition[2],
       ParticleVelocity[0], ParticleVelocity[1],
          ParticleVelocity[2],
       ParticleMass, ParticleAttribute[1], ParticleAttribute[0],
          ParticleAttribute[2], &RadiationData.IntegratedStarFormation);
 
  } // end: if (StarParticleFeedback == 2)
 
  if (StarParticleFeedback == 3) {

    //---- UNIGRID (NON-JEANS MASS) VERSION
 
      FORTRAN_NAME(star_feedback3)(
       GridDimension, GridDimension+1, GridDimension+2,
          BaryonField[DensNum], dmfield,
          BaryonField[TENum], BaryonField[GENum], BaryonField[Vel1Num],
          BaryonField[Vel2Num], BaryonField[Vel3Num], BaryonField[MetalNum],
          BaryonField[MetalNum+1], BaryonField[MetalNum+2],
       &DualEnergyFormalism, &MetallicityField, &MultiMetals, &HydroMethod,
       &dtFixed, BaryonField[NumberOfBaryonFields], &CellWidthTemp,
          &Time, &zred,
       &DensityUnits, &LengthUnits, &VelocityUnits, &TimeUnits,
          &StarEnergyToThermalFeedback, &StarMassEjectionFraction,
          &StarMetalYield,
       &NumberOfParticles,
          CellLeftEdge[0], CellLeftEdge[1], CellLeftEdge[2], &GhostZones,
       ParticlePosition[0], ParticlePosition[1],
          ParticlePosition[2],
       ParticleVelocity[0], ParticleVelocity[1],
          ParticleVelocity[2],
       ParticleMass, ParticleAttribute[1], ParticleAttribute[0],
          ParticleAttribute[2], &RadiationData.IntegratedStarFormation);
 
  } // end: if (StarParticleFeedback == 3)
 
  if (StarParticleFeedback == 4) {  

    //---- KRAVTSOV STAR FORMATION ALGORITHM

      FORTRAN_NAME(star_feedback4)(
       GridDimension, GridDimension+1, GridDimension+2,
          BaryonField[DensNum], dmfield, 
          BaryonField[TENum], BaryonField[GENum], BaryonField[Vel1Num],
          BaryonField[Vel2Num], BaryonField[Vel3Num], BaryonField[MetalNum],
       &DualEnergyFormalism, &MetallicityField, &HydroMethod, 
       &dtFixed, BaryonField[NumberOfBaryonFields], &CellWidthTemp, 
          &Time, &zred,
       &DensityUnits, &LengthUnits, &VelocityUnits, &TimeUnits,
          &StarEnergyToThermalFeedback, &StarMassEjectionFraction, 
          &StarMetalYield,
       &NumberOfParticles,
          CellLeftEdge[0], CellLeftEdge[1], CellLeftEdge[2], &GhostZones,
       ParticlePosition[0], ParticlePosition[1], 
          ParticlePosition[2], 
       ParticleVelocity[0], ParticleVelocity[1], 
          ParticleVelocity[2], 
       ParticleMass, ParticleAttribute[1], ParticleAttribute[0],
          ParticleAttribute[2], &RadiationData.IntegratedStarFormation);

  } // end: if (StarParticleFeedback == 4)


  if (StarParticleFeedback == 5) {  

    //---- SPRINGEL & HERNQUIST ALGORITHM

      FORTRAN_NAME(star_feedback5)(
       GridDimension, GridDimension+1, GridDimension+2,
          BaryonField[DensNum], dmfield, 
          BaryonField[TENum], BaryonField[GENum], BaryonField[Vel1Num],
          BaryonField[Vel2Num], BaryonField[Vel3Num], BaryonField[MetalNum],
       &DualEnergyFormalism, &MetallicityField, &HydroMethod, 
       &dtFixed, BaryonField[NumberOfBaryonFields], &CellWidthTemp, 
          &Time, &zred,
       &DensityUnits, &LengthUnits, &VelocityUnits, &TimeUnits,
          &StarEnergyToThermalFeedback, &StarMassEjectionFraction, 
          &StarMetalYield,
       &NumberOfParticles,
          CellLeftEdge[0], CellLeftEdge[1], CellLeftEdge[2], &GhostZones,
       ParticlePosition[0], ParticlePosition[1], 
          ParticlePosition[2], 
       ParticleVelocity[0], ParticleVelocity[1], 
          ParticleVelocity[2], 
       ParticleMass, ParticleAttribute[1], ParticleAttribute[0],
          ParticleAttribute[2], &RadiationData.IntegratedStarFormation);
      //fprintf(stderr,"After feedback is called");

  } // end: if (StarParticleFeedback == 5)

/* adding distributed starmaker and feedback */
  if (StarParticleFeedback == 10) {

    //---- UNIGRID (NON-JEANS MASS) VERSION
 
      FORTRAN_NAME(star_feedback10)(
       GridDimension, GridDimension+1, GridDimension+2,
          BaryonField[DensNum], dmfield,
          BaryonField[TENum], BaryonField[GENum], BaryonField[Vel1Num],
          BaryonField[Vel2Num], BaryonField[Vel3Num], BaryonField[MetalNum],
          BaryonField[MetalNum+1], BaryonField[MetalNum+2],
       &DualEnergyFormalism, &MetallicityField, &MultiMetals, &HydroMethod,
       &dtFixed, BaryonField[NumberOfBaryonFields], &CellWidthTemp,
          &Time, &zred,
       &DensityUnits, &LengthUnits, &VelocityUnits, &TimeUnits,
          &StarEnergyToThermalFeedback, &StarMassEjectionFraction,
          &StarMetalYield, &StarFeedbackDistRadius, &StarFeedbackDistCellStep, 
       &StarFeedbackDistTotalCells,
       &NumberOfParticles,
          CellLeftEdge[0], CellLeftEdge[1], CellLeftEdge[2], &GhostZones,
       ParticlePosition[0], ParticlePosition[1],
          ParticlePosition[2],
       ParticleVelocity[0], ParticleVelocity[1],
          ParticleVelocity[2],
       ParticleMass, ParticleAttribute[1], ParticleAttribute[0],
          ParticleAttribute[2], ParticleType, &RadiationData.IntegratedStarFormation);
 
  }


  /* Convert the species back from fractional densities to real densities. */
 
  for (field = 0; field < NumberOfBaryonFields; field++) {

    if (FieldType[field] >= ElectronDensity && 
#ifdef EMISSIVITY
        FieldType[field] <  GravPotential
#else
	FieldType[field] <  FieldUndefined
#endif

	) {
      for (k = GridStartIndex[2]; k <= GridEndIndex[2]; k++) {
	for (j = GridStartIndex[1]; j <= GridEndIndex[1]; j++) {
	  index = (k*GridDimension[1] + j)*GridDimension[0] +
	    GridStartIndex[0];
	  for (i = GridStartIndex[0]; i <= GridEndIndex[0]; i++, index++) {
	    BaryonField[field][index] *= BaryonField[DensNum][index];
	  }
	}
      }

    }

  }
 
  /* Clean up. */
 
  delete [] temperature;
  delete [] dmfield;
  delete [] BaryonField[NumberOfBaryonFields];   // refinement flag field
  BaryonField[NumberOfBaryonFields] = NULL;
 
  if (debug) printf("StarParticle: end\n");

  return SUCCESS;
}
