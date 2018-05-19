/***********************************************************************
/
/  SETS THE DEFAULT GLOBAL VALUES
/
/  written by: Greg Bryan
/  date:       November, 1994
/  modified1:  Robert Harkness
/  date:       February 29th, 2008
/  modified2:  S.W. Skillman
/  date:       January, 2009
/  modified3:  Robert Harkness
/  date:       January, 2010
/              Set DensityType and Chaining Mesh parameters
/
/  PURPOSE:
/
/  RETURNS: SUCCESS or FAIL
/
************************************************************************/
 
// This routine intializes a new simulation based on the parameter file.
//
 
#include <string.h>
#include <stdio.h>
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "TopGridData.h"
#include "StarParticleData.h"

int SetFieldDensityType(void);
int PrintFieldDensityType(void);
 
/* character strings */
 
char DefaultDimUnits[] = "cm";
char *DefaultDimLabel[] = {"x", "y", "z"};
 
char DefaultRestartName[] = "restart";
char DefaultDataName[] = "data";
char DefaultHistoryName[] = "history";
char DefaultRedshiftName[] = "RedshiftOutput";
char DefaultMovieName[] = "MovieOutput";
char DefaultNewMovieName[] = "MoviePack";
char DefaultTracerParticleName[] = "TracerOutput";
 
char DefaultRestartDir[] = "RS";
char DefaultDataDir[] = "DD";
char DefaultHistoryDir[] = "HD";
char DefaultRedshiftDir[] = "RD";
char DefaultMovieDir[] = "MD";
char DefaultTracerParticleDir[] = "TD";
 
 
 
 
int SetDefaultGlobalValues(TopGridData &MetaData)
{
 
  /* declarations */
 
  const float Pi = 3.14159;
  int dim, i;
 
  /* set the default MetaData values. */
 
  MetaData.CycleNumber     = 0;
  MetaData.SubcycleNumber  = 0;
  MetaData.Time            = 0.0;
  MetaData.CPUTime         = 0.0;
 
  MetaData.StopTime        = FLOAT_UNDEFINED;  // This must be set be the user
  MetaData.StopCycle       = 10000;            // 10000 timesteps
  MetaData.StopSteps       = 10000;            // 10000 timesteps
  MetaData.StopCPUTime     = 100.0*3600.0;     // 100 hours
 
  MetaData.TimeLastRestartDump = 0.0;
  MetaData.dtRestartDump       = 5.0*3600.0;   // every 5 hours
  MetaData.TimeLastDataDump    = FLOAT_UNDEFINED;
  MetaData.dtDataDump          = 0.0;
  MetaData.TimeLastHistoryDump = FLOAT_UNDEFINED;
  MetaData.dtHistoryDump       = 0.0;
  MetaData.TimeLastMovieDump   = FLOAT_UNDEFINED;
  MetaData.dtMovieDump         = 0.0;
  MetaData.TimeLastTracerParticleDump = FLOAT_UNDEFINED;
  MetaData.dtTracerParticleDump       = 0.0;
 
  MetaData.CycleLastRestartDump = 0;
  MetaData.CycleSkipRestartDump = 0;
  MetaData.CycleLastDataDump    = INT_UNDEFINED;
  MetaData.CycleSkipDataDump    = 0;
  MetaData.SubcycleLastDataDump = INT_UNDEFINED;
  MetaData.SubcycleSkipDataDump = 0;
  MetaData.CycleLastHistoryDump = INT_UNDEFINED;
  MetaData.CycleSkipHistoryDump = 0;
  MetaData.CycleSkipGlobalDataDump = 0; //AK
 
  MetaData.OutputFirstTimeAtLevel = 0; // zero is off
  MetaData.StopFirstTimeAtLevel   = 0; // zero is off
 
  MetaData.RestartDumpNumber   = 0;            // starting restart id number
  MetaData.RestartDumpName     = DefaultRestartName;
  MetaData.RestartDumpDir      = DefaultRestartDir;
  MetaData.DataDumpNumber      = 0;
  MetaData.DataDumpName        = DefaultDataName;
  MetaData.DataDumpDir         = DefaultDataDir;
  MetaData.HistoryDumpNumber   = 0;
  MetaData.HistoryDumpName     = DefaultHistoryName;
  MetaData.HistoryDumpDir      = DefaultHistoryDir;
  MetaData.MovieDumpNumber     = 0;
  MetaData.MovieDumpName       = DefaultMovieName;
  MetaData.MovieDumpDir        = DefaultMovieDir;
  MetaData.TracerParticleDumpNumber = 0;
  MetaData.TracerParticleDumpName   = DefaultTracerParticleName;
  MetaData.TracerParticleDumpDir    = DefaultTracerParticleDir;
//MetaData.RedshiftDumpNumber  = 0;
  MetaData.RedshiftDumpName    = DefaultRedshiftName;
  MetaData.RedshiftDumpDir     = DefaultRedshiftDir;
 
  MetaData.LocalDir            = NULL;
  MetaData.GlobalDir           = NULL;
 
  for (i = 0; i < MAX_TIME_ACTIONS; i++) {
    TimeActionType[i]      = 0;
    TimeActionRedshift[i]  = -1;
    TimeActionTime[i]      = 0;
    TimeActionParameter[i] = FLOAT_UNDEFINED;
  }
 
  for (i = 0; i < MAX_CUBE_DUMPS; i++) {
    CubeDumps[i] = NULL;
  }
 
  MetaData.StaticHierarchy     = TRUE;
 
  MetaData.TopGridRank = INT_UNDEFINED;
  for (dim = 0; dim < MAX_DIMENSION; dim++) {
    MetaData.TopGridDims[dim]                = INT_UNDEFINED;
    MetaData.LeftFaceBoundaryCondition[dim]  = reflecting;
    MetaData.RightFaceBoundaryCondition[dim] = reflecting;
  }
  MetaData.BoundaryConditionName = NULL;
 
  MetaData.GravityBoundary        = TopGridPeriodic;

#ifdef RAD_HYDRO
  MetaData.RadHydroParameterFname = NULL;
#endif
 
  MetaData.ParticleBoundaryType   = periodic;  // only one implemented!
  MetaData.NumberOfParticles      = 0;         // no particles
 
  MetaData.CourantSafetyNumber    = 0.6;
  MetaData.PPMFlatteningParameter = 0;    // off
  MetaData.PPMDiffusionParameter  = 0;    // off
  MetaData.PPMSteepeningParameter = 0;    // off
 
  /* set the default global data. */
                                                 // Debug flag set in main
  ProblemType               = 0;                 // None
  HydroMethod               = PPM_DirectEuler;   //
  huge_number               = 1.0e+20;
  tiny_number               = 1.0e-20;
  Gamma                     = 5.0/3.0;           // 5/3
  PressureFree              = FALSE;             // use pressure (duh)
  RefineBy                  = 4;                 // Refinement factor
  MaximumRefinementLevel    = 2;                 // three levels (w/ topgrid)
  MaximumGravityRefinementLevel = INT_UNDEFINED;
  MaximumParticleRefinementLevel = -1;            // unused if negative
  MustRefineRegionMinRefinementLevel = -1;        // unused if negative
  MetallicityRefinementMinLevel = -1;
  MetallicityRefinementMinMetallicity = 1.0e-5;
  FluxCorrection            = TRUE;
  InterpolationMethod       = SecondOrderA;      // ?
  ConservativeInterpolation = TRUE;              // true for ppm
  MinimumEfficiency         = 0.2;               // between 0-1, usually ~0.1
  MinimumSubgridEdge        = 4;                 // min for acceptable subgrid
  MaximumSubgridSize        = 2000;              // max for acceptable subgrid
  NumberOfBufferZones       = 1;
 
  for (i = 0; i < MAX_FLAGGING_METHODS; i++) {
    CellFlaggingMethod[i]       = INT_UNDEFINED;
    MinimumMassForRefinement[i] = FLOAT_UNDEFINED;   // usually set by:
    MinimumOverDensityForRefinement[i]       = 1.5;
    MinimumMassForRefinementLevelExponent[i] = 0;
  }
 
  for (dim = 0; dim < MAX_DIMENSION; dim++) {
    DomainLeftEdge[dim]             = 0.0;
    DomainRightEdge[dim]            = 1.0;
    GridVelocity[dim]               = 0.0;
    DimLabels[dim]                  = DefaultDimLabel[dim];
    DimUnits[dim]                   = DefaultDimUnits;
    RefineRegionLeftEdge[dim]       = FLOAT_UNDEFINED;
    RefineRegionRightEdge[dim]      = FLOAT_UNDEFINED;
    MetaData.MovieRegionLeftEdge[dim]  = FLOAT_UNDEFINED;
    MetaData.MovieRegionRightEdge[dim] = FLOAT_UNDEFINED;
    MetaData.NewMovieLeftEdge[dim]  = 0.0;
    MetaData.NewMovieRightEdge[dim] = 1.0;
    PointSourceGravityPosition[dim] = 0.0;
    MustRefineRegionLeftEdge[dim] = 0.0;
    MustRefineRegionRightEdge[dim] = 1.0;
  }
 
  for (i = 0; i < MAX_STATIC_REGIONS; i++)
    StaticRefineRegionLevel[i] = INT_UNDEFINED;
 
  ParallelRootGridIO          = FALSE;
  ParallelParticleIO          = FALSE;
  Unigrid                     = FALSE;
  PartitionNestedGrids        = FALSE;
  ExtractFieldsOnly           = TRUE;
  SRBprefix                   = NULL;
  SRBcwd                      = new char[MAX_LINE_LENGTH];
  PackedAMR                   = TRUE;

  ExternalBoundaryIO          = FALSE;
  ExternalBoundaryTypeIO      = FALSE;
  ExternalBoundaryValueIO     = FALSE;
  SimpleConstantBoundary      = FALSE;

  StaticChainingMeshDimension = 16;
  DynamicChainingMeshDimension = 16;

  LoadBalance0                = 0;
  LoadBalance1                = 0;
  LoadBalance2                = 0;

  debug1                      = 0;
  debug2                      = 0;

  TracerParticleOn            = 0;
 
  CubeDumpEnabled             = 0;

#ifdef STAGE_INPUT
  StageInput                  = 0;
#endif

  First_Pass                  = 0;

  MemoryLimit                 = 4000000000;
 
  UniformGravity              = FALSE;             // off
  UniformGravityDirection     = 0;                 // x-direction
  UniformGravityConstant      = 1.0;
 
  PointSourceGravity           = FALSE;             // off
  PointSourceGravityConstant   = 1.0;
  PointSourceGravityCoreRadius = 0.0;

  SelfGravity                 = FALSE;             // off
  CopyGravPotential           = FALSE;             // off
  PotentialIterations         = 0;                 // ~4 is reasonable
  GravitationalConstant       = 4*Pi;              // G = 1
  S2ParticleSize              = 3.0;               // ~3 is reasonable
  GravityResolution           = 1.0;               // equivalent to grid
  ComputePotential            = FALSE;
  WritePotential              = FALSE;
  BaryonSelfGravityApproximation = TRUE;           // less accurate but faster
 
  GreensFunctionMaxNumber     = 1;                 // only one at a time
  GreensFunctionMaxSize       = 1;                 // not used yet
 
  DualEnergyFormalism         = FALSE;             // off
  DualEnergyFormalismEta1     = 0.001;             // typical 0.001
  DualEnergyFormalismEta2     = 0.1;               // 0.08-0.1
  ParticleCourantSafetyNumber = 0.5;
  RootGridCourantSafetyNumber = 1.0;
  RandomForcing               = FALSE;             // off //AK
  RandomForcingEdot           = -1.0;              //AK
  RandomForcingMachNumber     = 0.0;               //AK
  RadiativeCooling            = FALSE;             // off
  RadiativeTransfer           = FALSE;             // off
  RadiativeTransferHydrogenOnly = FALSE;           // off
#ifdef RAD_HYDRO
  ImplicitProblem             = 0;                 // off
  RadiationHydrodynamics      = 0;                 // off
  StarMakerEmissivityField    = 0;                 // off
#endif
  UV_parameter                = 0.0;               // Razoumov and Norman Ap.J. 2002 June 20 5.0e-04
  Xr_parameter                = 0.0;               
  GadgetEquilibriumCooling    = FALSE;             // off
  MultiSpecies                = FALSE;             // off
  CRModel                     = 0;                 // off
  ShockMethod                 = 0;                 // temperature unsplit
  ShockTemperatureFloor       = 1.0;               // Set to 1K
  StorePreShockFields         = 0;
  GloverChemistryModel        = 0;                 // 0ff
  GloverRadiationBackground   = 0;
  GloverOpticalDepth          = 0;
  RadiationFieldType          = 0;
  RadiationFieldLevelRecompute = 0;
  AdjustUVBackground          = 1;
  SetUVBAmplitude             = 1.0;
  SetHeIIHeatingScale         = 1.8;
  CoolData.alpha0             = 1.5;               // radiation spectral slope
  CoolData.f3                 = 1.0e-21;           // radiation normalization
  CoolData.ParameterFilename  = NULL;

  MetalCooling                = 0;                 // off
  CloudyCoolingData.CloudyCoolingGridRank = 0;
  CloudyCoolingData.CloudyCoolingGridFile = "";
  CloudyCoolingData.IncludeCloudyHeating = 0;
  CloudyCoolingData.IncludeCloudyMMW = 0;
  CloudyCoolingData.CMBTemperatureFloor         = 1; // use CMB floor.
  CloudyCoolingData.CloudyMetallicityNormalization = 0.018477;      // calculated using Cloudy 07.02 abundances
  CloudyCoolingData.CloudyElectronFractionFactor = 9.153959e-3;     // calculated using Cloudy 07.02 abundances

  OutputCoolingTime = 0;
  OutputTemperature = 0;

  ZEUSLinearArtificialViscosity    = 0.0;
  ZEUSQuadraticArtificialViscosity = 2.0;
  UseMinimumPressureSupport        = FALSE;
  MinimumPressureSupportParameter  = 100.0;
 
  MinimumSlopeForRefinement        = 0.3;          // 30% change in value
  MinimumShearForRefinement        = 1.0;          //AK
  MinimumPressureJumpForRefinement = 0.33;         // As in PPM method paper
  MinimumEnergyRatioForRefinement  = 0.1;          // conservative!
  RefineByJeansLengthSafetyFactor  = 4.0;
  MustRefineParticlesRefineToLevel = 0;
  ComovingCoordinates              = FALSE;        // No comoving coordinates
  StarParticleCreation             = FALSE;
  StarParticleFeedback             = FALSE;
  StarMakerOverDensityThreshold    = 100;          // times mean total density
  StarMakerMassEfficiency          = 1;
  StarMakerMinimumMass             = 1.0e9;        // in solar masses
  StarMakerMinimumDynamicalTime    = 1.0e6;        // in years
  StarMassEjectionFraction         = 0.25;
  StarMetalYield                   = 0.02;
  StarEnergyToThermalFeedback      = 1.0e-5;
  StarEnergyToStellarUV            = 3.0e-6;
  StarEnergyToQuasarUV             = 5.0e-6;
/* adding distributed starmaker and feedback */
  StarFeedbackDistRadius           = 0;
  StarFeedbackDistCellStep         = 0;
  StarFeedbackDistTotalCells       = 1;
  MultiMetals                      = FALSE;
  NumberOfParticleAttributes       = INT_UNDEFINED;

  ParticleTypeInFile  = FALSE;
 
  for (i = 0; i<MAX_MOVIE_FIELDS; i++)
    MovieDataField[i] = INT_UNDEFINED;
  MovieSkipTimestep = INT_UNDEFINED;
  NewMovieName = DefaultNewMovieName;
  NewMovieDumpNumber = 0;
  NewMovieEntries = 0;
  MovieEntriesPP = new long[NumberOfProcessors];
  MaxMovieFilenum = 0;
  for (i = 0; i<NumberOfProcessors; i++)
    MovieEntriesPP[i] = 0;
  NewMovieParticleOn = FALSE;

  ran1_init = 0;

  RiemannSolver              = TwoShock; 
  RiemannSolverFallback      = 0;
  ReconstructionMethod       = PPM;
  PositiveReconstruction     = 0;
  ConservativeReconstruction = 0;

  /* test problem values */
  TestProblemData.HydrogenFractionByMass = 0.76;

  /* The DToHRatio is by mass in the code, so multiply by 2. */
  TestProblemData.DeuteriumToHydrogenRatio = 2.0*3.4e-5; // Burles & Tytler 1998

  // multispecies default values assume completely neutral gas with primordial D/H ratio
  TestProblemData.MultiSpecies = 0;
  TestProblemData.HI_Fraction = 1.0;
  TestProblemData.HII_Fraction = tiny_number;
  TestProblemData.HeI_Fraction = 1.0;
  TestProblemData.HeII_Fraction = tiny_number;
  TestProblemData.HeIII_Fraction = tiny_number;
  TestProblemData.HM_Fraction = tiny_number;
  TestProblemData.H2I_Fraction = tiny_number;
  TestProblemData.H2II_Fraction = tiny_number;
  TestProblemData.DI_Fraction = 2.0*3.4e-5;
  TestProblemData.DII_Fraction = tiny_number;
  TestProblemData.HDI_Fraction = tiny_number;

  // This is for ionized gas (within a supernova blast, for example)
  TestProblemData.HI_Fraction_Inner = 1.0;
  TestProblemData.HII_Fraction_Inner = tiny_number;
  TestProblemData.HeI_Fraction_Inner = 1.0;
  TestProblemData.HeII_Fraction_Inner = tiny_number;
  TestProblemData.HeIII_Fraction_Inner = tiny_number;
  TestProblemData.HM_Fraction_Inner = tiny_number;
  TestProblemData.H2I_Fraction_Inner = tiny_number;
  TestProblemData.H2II_Fraction_Inner = tiny_number;
  TestProblemData.DI_Fraction_Inner = 2.0*3.4e-5;
  TestProblemData.DII_Fraction_Inner = tiny_number;
  TestProblemData.HDI_Fraction_Inner = tiny_number;

  TestProblemData.UseMetallicityField = 0;
  TestProblemData.MetallicityField_Fraction = tiny_number;

  TestProblemData.UseMassInjection = 0;
  TestProblemData.InitialHydrogenMass = tiny_number;
  TestProblemData.InitialDeuteriumMass = tiny_number;
  TestProblemData.InitialHeliumMass = tiny_number;
  TestProblemData.InitialMetalMass = tiny_number;

  TestProblemData.MultiMetals = 0;
  TestProblemData.MultiMetalsField1_Fraction = tiny_number;
  TestProblemData.MultiMetalsField2_Fraction = tiny_number;

  TestProblemData.GloverChemistryModel = 0;
  // This is for the gas in the surrounding medium, for the blast wave problem.
  TestProblemData.CI_Fraction = tiny_number;
  TestProblemData.CII_Fraction = tiny_number;
  TestProblemData.OI_Fraction = tiny_number;
  TestProblemData.OII_Fraction = tiny_number;
  TestProblemData.SiI_Fraction = tiny_number;
  TestProblemData.SiII_Fraction = tiny_number;
  TestProblemData.SiIII_Fraction = tiny_number;
  TestProblemData.CHI_Fraction = tiny_number;
  TestProblemData.CH2I_Fraction = tiny_number;
  TestProblemData.CH3II_Fraction = tiny_number;
  TestProblemData.C2I_Fraction = tiny_number;
  TestProblemData.COI_Fraction = tiny_number;
  TestProblemData.HCOII_Fraction = tiny_number;
  TestProblemData.OHI_Fraction = tiny_number;
  TestProblemData.H2OI_Fraction = tiny_number;
  TestProblemData.O2I_Fraction = tiny_number;

  // This is for the gas in the region where the blast wave energy is injected
  TestProblemData.CI_Fraction_Inner = tiny_number;
  TestProblemData.CII_Fraction_Inner = tiny_number;
  TestProblemData.OI_Fraction_Inner = tiny_number;
  TestProblemData.OII_Fraction_Inner = tiny_number;
  TestProblemData.SiI_Fraction_Inner = tiny_number;
  TestProblemData.SiII_Fraction_Inner = tiny_number;
  TestProblemData.SiIII_Fraction_Inner = tiny_number;
  TestProblemData.CHI_Fraction_Inner = tiny_number;
  TestProblemData.CH2I_Fraction_Inner = tiny_number;
  TestProblemData.CH3II_Fraction_Inner = tiny_number;
  TestProblemData.C2I_Fraction_Inner = tiny_number;
  TestProblemData.COI_Fraction_Inner = tiny_number;
  TestProblemData.HCOII_Fraction_Inner = tiny_number;
  TestProblemData.OHI_Fraction_Inner = tiny_number;
  TestProblemData.H2OI_Fraction_Inner = tiny_number;
  TestProblemData.O2I_Fraction_Inner = tiny_number;

  TestProblemData.MinimumHNumberDensity = 1;
  TestProblemData.MaximumHNumberDensity = 1e6;
  TestProblemData.MinimumMetallicity    = 1e-6;
  TestProblemData.MaximumMetallicity    = 1;
  TestProblemData.MinimumTemperature    = 10;
  TestProblemData.MaximumTemperature    = 1e7;
  TestProblemData.ResetEnergies         = 1;

  // This should only be false for analysis.
  // It could also be used (cautiously) for other purposes.
  LoadGridDataAtStart = TRUE;

  IsothermalSoundSpeed = 1.0;
  RefineByJeansLengthUnits = 0;

  SetFieldDensityType();
  PrintFieldDensityType();

  return SUCCESS;
}
