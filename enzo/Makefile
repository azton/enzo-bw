#=======================================================================
#  MAKEFILE FOR ENZO-R
#  Robert Harkness
#=======================================================================

VERSION := "03 January 2011"
DATE := $(shell date)
CWD := $(shell pwd)

# Settings - change at will
ENZO_SYSTEM := reynolds
ENZO_COMPILER := intel
#ENZO_SYSTEM := gso_macports
#ENZO_COMPILER := gnu
#ENZO_SYSTEM := bluedrop
#ENZO_COMPILER := ibmxl
#ENZO_SYSTEM := kraken
#ENZO_COMPILER := pgi

#  MODE is hybrid or mpi
ENZO_MODE := mpi
ENZO_SPEC := _test
ENZO_AMR := no
ENZO_RHD := yes
ENZO_TAG := LARGE_TAGS
ENZO_INT := LARGE_INTS
HYDRO := WISE2

MAX_GRIDS := 3000
MAX_BARYONS := 25
MAX_COLORS := 25
MAX_CUBE_PLUS_7 = 1031
MAX_CPU_PER_NODE := 2

# Don't change anything below this if you don't have to!

# Precision of initial conditions
INIT_PPR := INITS64

# Optimizations - don't change these unless you know why
OOC_BDRY := yes
UXP_OPT := UNIGRID_TRANSPOSE
FFT_OPT := USE_FFTE_FFT
IBM_TASKMAP := no

# Legacy Options
ENZO_SRB := no
ENZO_END := HDF5_BE
ISO_GRAV := NO_ISO_GRAV
JB_PERF := NO_USE_JBPERF

# Special Options

# AMR_OPT := -D SIB1 -D SIB2 -D SIB3 -D SIB4 -D SIB5 -D STATIC_SIBLING_LIST -D FLUX_FIX -D SAB \
#            -D JB_OPT_FAST_NEIGHBOR_SEARCH -D USE_DT_LIMIT -D FORCE_BUFFER_PURGE -D FORCE_MSG_PROGRESS -D ENABLE_LOAD_BALANCE
#
# These settings used with the 500 Mpc 512^3 AMR, i.e. disable new flux correction

ifeq ($(ENZO_RHD),yes)
  RHD_OPT := -DRAD_HYDRO -DEMISSIVITY -DUSE_HYPRE
else
  RHD_OPT :=
endif

ifeq ($(ENZO_AMR),yes)
  PACK_AMR := yes
  PACK_MEM := yes
  ENZO_XPR := r8
  ENZO_PPR := p8
  PREC_IO := yes
  IO_PR := IO_64
  AMR_OPT := -D SIB1 -D SIB2 -D SIB3 -D SIB4 -D SIB5 -D NO_FLUX_FIX -D NO_SAB -D NO_DC_OPT_SIBSUB -D JB_OPT_FAST_NEIGHBOR_SEARCH -D DONT_USE_DT_LIMIT -D FORCE_BUFFER_PURGE -D FORCE_MSG_PROGRESS -D DISABLE_LOAD_BALANCE -D MISCOUNT -D NO_RATE_AND_COOL
else
  PACK_AMR := no
  PACK_MEM := no
  ifeq ($(ENZO_RHD),yes)
    ENZO_XPR := r8
    ENZO_PPR := p8
    PREC_IO := yes
    IO_PR := IO_64
  else
    ENZO_XPR := r8
    ENZO_PPR := p8
    PREC_IO := no
    IO_PR := IO_32
  endif
  AMR_OPT := -D SIB1 -D SIB2 -D SIB3 -D SIB4 -D SIB5 -D STATIC_SIBLING_LIST -D NO_FLUX_FIX -D NO_SAB -D NO_DC_OPT_SIBSUB -D JB_OPT_FAST_NEIGHBOR_SEARCH -D DONT_USE_DT_LIMIT -D FORCE_BUFFER_PURGE -D FORCE_MSG_PROGRESS -D DISABLE_LOAD_BALANCE -D MISCOUNT -D NO_RATE_AND_COOL
endif


ifeq ($(PACK_AMR),yes)
  IO_PACK := USE_HDF5_GROUPS
  IO_OPT := SINGLE_HDF5_OPEN_ON_INPUT
  ifeq ($(PACK_MEM),yes)
    IO_MEM_BUFF_OUT := USE_HDF5_OUTPUT_BUFFERING
    IO_MEM_BUFF_INP := DONT_USE_HDF5_INPUT_BUFFERING
  else
    IO_MEM_BUFF_OUT := NO_HDF5_OUTPUT_BUFFERING
    IO_MEM_BUFF_INP := NO_HDF5_INPUT_BUFFERING
  endif
else
  IO_PACK := NO_HDF5_GROUPS
  IO_OPT := NO_SINGLE_HDF5_OPEN_ON_INPUT
  IO_MEM_BUFF_OUT := NO_HDF5_OUTPUT_BUFFERING
  IO_MEM_BUFF_INP := NO_HDF5_INPUT_BUFFERING
endif

ifeq ($(IBM_TASKMAP),yes)
  TASK_MAP := TASKMAP
  USE_LOADL := DONT_USE_LL
else
  TASK_MAP := NO_TASKMAP
  USE_LOADL := DONT_USE_LL
endif

ifeq ($(OOC_BDRY),yes)
  OOC_OPT := OOC_BOUNDARY
else
  OOC_OPT := NO_OOC_BOUNDARY
endif

#-------------------  object files  ------------------------------------


F90X = 	select_fft.o					\
	fft90.o

XOBJ =	AdiabaticExpansionInitialize.o			\
	SetFieldColorType.o				\
	PrintFieldColorType.o				\
	SetFieldDensityType.o				\
	PrintFieldDensityType.o				\
	AssignGridToTaskMap.o				\
	CheckEnergyConservation.o			\
	CheckForCycleDump.o				\
	CheckForDump.o					\
	CheckForHalt.o					\
	CheckForOutput.o				\
	CheckForTimeAction.o				\
	CollapseTestInitialize.o			\
	CommunicationBroadcastValue.o			\
	CommunicationBufferedSend.o			\
	CommunicationCombineGrids.o			\
	CommunicationInitialize.o			\
	CommunicationLoadBalanceGrids.o			\
	CommunicationParallelFFT.o			\
	CommunicationPartitionGrid.o			\
	CommunicationReceiveFluxes.o			\
	CommunicationSendFluxes.o			\
	CommunicationShareGrids.o			\
	CommunicationTransferParticles.o		\
	CommunicationTranspose.o			\
	CommunicationUpdateStarParticleCount.o		\
	CommunicationUtilities.o			\
	ComputePotentialFieldLevelZero.o		\
	ComputeTable.o					\
	ComputeRandomForcingNormalization.o		\
	Continue.o					\
	CoolingTestInitialize.o				\
	CopyOverlappingParticleMassFields.o		\
	CopyOverlappingZones.o				\
	CosmologyComputeExpansionFactor.o		\
	CosmologyComputeExpansionTimestep.o		\
	CosmologyComputeTimeFromRedshift.o		\
	CosmologyGetUnits.o				\
	CosmologyReadParameters.o			\
	CosmologySimulationInitialize.o			\
	NestedCosmologySimulationInitialize.o		\
	CosmologyWriteParameters.o			\
	CreateSUBlingList.o				\
	DeleteSUBlingList.o				\
	DeleteFluxes.o					\
	DepositBaryons.o				\
	DepositParticleMassField.o			\
	DoubleMachInitialize.o				\
	Enzo_Dims_create.o				\
	Mpich_V1_Dims_create.o				\
	EvolveHierarchy.o				\
	ExternalBoundary_AppendForcingToBaryonFields.o	\
	ExternalBoundary_DetachForcingFromBaryonFields.o	\
	ExternalBoundary_IdentifyPhysicalQuantities.o	\
	ExternalBoundary_InitializeExternalBoundaryFaceIO.o	\
	ExternalBoundary_Prepare.o			\
	ExternalBoundary_SetDoubleMachBoundary.o	\
	ExternalBoundary_SetExternalBoundaryIO.o	\
	ExternalBoundary_SetExternalBoundaryParticles.o	\
	ExternalBoundary_SetShockPoolBoundary.o		\
	ExternalBoundary_SetWavePoolBoundary.o		\
	ExternalBoundary_constructor.o			\
	FastFourierTransform.o				\
	FastFourierTransformPrepareComplex.o		\
	FastFourierTransformSGIMATH.o			\
	FastSiblingLocatorInitialize.o			\
	FastSiblingLocatorInitializeStaticChainingMesh.o \
	FastSiblingLocatorFinalize.o			\
	FindCube.o					\
	FindField.o					\
	FindSubgrids.o					\
	FreezeRateData.o				\
	Gadget.o					\
	GadgetInitCoolMemory.o				\
	GadgetIonizeParamsTable.o			\
	GadgetMakeCoolingTable.o			\
	GadgetReadIonizeParams.o			\
	GalaxySimulationInitialize.o			\
	GetNodeFreeMemory.o				\
	GetUnits.o					\
	GravityBdryExchange.o				\
	GravityEquilibriumTestInitialize.o		\
	Grid_AccelerationBoundaryRoutines.o		\
	Grid_AddBaryonsToGravitatingMassField.o		\
	Grid_AddExternalAcceleration.o			\
	Grid_AddFieldMassToMassFlaggingField.o		\
	Grid_AddOverlappingParticleMassField.o		\
	Grid_AddRandomForcing.o				\
	Grid_AddToBoundaryFluxes.o			\
	Grid_AllocateGrids.o				\
	Grid_AnalyzeTrackPeaks.o			\
	Grid_AppendForcingToBaryonFields.o		\
	Grid_ApplyTimeAction.o				\
	Grid_CalculateAngularMomentum.o			\
	Grid_CheckForExternalReflections.o		\
	Grid_CheckForOverlap.o				\
	Grid_CheckForPossibleOverlap.o			\
	Grid_CheckForSharedFace.o			\
	Grid_CleanUp.o					\
	Grid_CleanUpMovedParticles.o			\
	Grid_ClearBoundaryFluxes.o			\
	Grid_ClearFlaggingField.o			\
	Grid_ClearGravitatingMassField.o		\
	Grid_ClearGravitatingMassFieldParticles.o	\
	Grid_ClearMassFlaggingField.o			\
	Grid_ClearParticleAccelerations.o		\
	Grid_CollapseTestInitializeGrid.o		\
	Grid_CollectGridInformation.o			\
	Grid_CommunicationMoveGrid.o			\
	Grid_IdentifyRadiativeTransferFields.o		\
	Grid_Irecv_CommunicationSendParticles.o		\
	Grid_Irecv_CommunicationSendRegion.o		\
	Grid_Irecv_CommunicationReceiveRegion.o		\
	Grid_CommunicationTransferParticles.o		\
	Grid_ComovingExpansionTerms.o			\
	Grid_ComovingGravitySourceTerm.o		\
	Grid_ComputeAccelerationField.o			\
	Grid_ComputeAccelerationFieldExternal.o		\
	Grid_ComputeAccelerations.o			\
	Grid_ComputeElementalDensity.o			\
	Grid_ComputeEnergy.o				\
	Grid_ComputeExternalNohBoundary.o		\
	Grid_ComputeGammaField.o			\
	Grid_ComputePressure.o				\
	Grid_ComputePressureDualEnergyFormalism.o	\
	Grid_ComputeTemperatureField.o			\
	Grid_ComputeTimeStep.o				\
	Grid_ComputeXrayEmissivity.o			\
	Grid_ConvertTotalEnergyToGasEnergy.o		\
	Grid_CoolingTestInitializeGrid.o		\
	Grid_CopyOverlappingMassField.o			\
	Grid_CopyParentToGravitatingFieldBoundary.o	\
	Grid_CopyPotentialField.o			\
	Grid_CopyPotentialToBaryonField.o		\
	Grid_CopyZonesFromGrid.o			\
	Grid_CopyZonesFromGridCountOnly.o		\
	Grid_CountGridsOnLevel.o			\
	Grid_DebugCheck.o				\
	Grid_DeleteAllFields.o				\
	Grid_DeleteFlaggingField.o			\
	Grid_DetachForcingFromBaryonFields.o		\
	Grid_Irecv_DepositBaryons.o			\
	Grid_DepositParticlePositions.o			\
	Grid_DepositPositions.o				\
	Grid_DoubleMachInitializeGrid.o			\
	Grid_FastSiblingLocatorAddGrid.o		\
	Grid_FastSiblingLocatorFindSiblings.o		\
	Grid_FindShocks.o				\
	Grid_FinishFFT.o				\
	Grid_FlagBufferZones.o				\
	Grid_FlagCellsToBeRefinedByCoolingTime.o	\
	Grid_FlagCellsToBeRefinedByJeansLength.o	\
	Grid_FlagCellsToBeRefinedByMass.o		\
	Grid_FlagCellsToBeRefinedByMetallicity.o	\
	Grid_FlagCellsToBeRefinedByMustRefineParticles.o	\
	Grid_FlagCellsToBeRefinedByMustRefineRegion.o	\
	Grid_FlagCellsToBeRefinedBySlope.o		\
	Grid_FlagCellsToBeRefinedByShear.o		\
	Grid_FlagCellsToBeRefinedByShocks.o		\
	Grid_GadgetCalculateCooling.o			\
	Grid_GadgetComputeTemperatureDEF.o		\
	Grid_GadgetComputeTemperature.o			\
	Grid_Gadgetconvertutotemp.o			\
	Grid_GadgetCoolingRateFromU.o			\
	Grid_GadgetCoolingRate.o			\
	Grid_GadgetCoolingTime.o			\
	Grid_GadgetEquilibriumCooling.o			\
	Grid_Gadgetfindabundancesandrates.o		\
	Grid_GalaxySimulationInitializeGrid.o		\
	Grid_GravityEquilibriumTestInitializeGrid.o	\
	Grid_IdentifyCRSpeciesFields.o			\
	Grid_IdentifyGloverSpeciesFields.o		\
	Grid_IdentifyNewSubgrids.o			\
	Grid_IdentifyNewSubgridsSmall.o			\
	Grid_IdentifyPhysicalQuantities.o		\
	Grid_IdentifySpeciesFields.o			\
	Grid_ImplosionInitializeGrid.o			\
	Grid_InheritProperties.o			\
	Grid_InitializeGravitatingMassField.o		\
	Grid_InitializeGravitatingMassFieldParticles.o	\
	Grid_InitializeUniformGrid.o			\
	Grid_InterpolateAccelerations.o			\
	Grid_InterpolateBoundaryFromParent.o		\
	Grid_InterpolateFieldValues.o			\
	Grid_InterpolateParticlePositions.o		\
	Grid_InterpolatePositions.o			\
	Grid_KHInitializeGrid.o				\
	Grid_MoveAllParticles.o				\
	Grid_MoveSubgridParticles.o			\
	Grid_MoveSubgridParticlesFast.o			\
	Grid_NohInitializeGrid.o			\
	Grid_OutputAsParticleData.o			\
	Grid_OutputGridMovieData.o			\
	Grid_OutputStarParticleInformation.o		\
	Grid_PrepareBoundaryFluxes.o			\
	Grid_PrepareFFT.o				\
	Grid_PrepareGreensFunction.o			\
	Grid_PreparePeriodicGreensFunction.o		\
	Grid_PrepareGrid.o				\
	Grid_PrepareGridDerivedQuantities.o		\
	Grid_PreparePotentialField.o			\
	Grid_PrepareRandomForcingNormalization.o	\
	Grid_PressurelessCollapseInitialize.o		\
	Grid_ProjectToPlane.o				\
	Grid_ProtostellarCollapseInitializeGrid.o	\
	Grid_RadiationComputeDensities.o		\
	Grid_RadiatingShockInitializeGrid.o		\
	Grid_ReturnFluxDims.o				\
	Grid_ReturnGridInfo.o				\
	Grid_RemoveForcingFromBaryonFields.o		\
	Grid_RestoreEnergyConsistency.o			\
	Grid_RotatingCylinderInitialize.o		\
	Grid_SedovBlastInitializeGrid.o			\
	Grid_SedovBlastInitializeGrid3D.o		\
	Grid_SetExternalBoundaryValues.o		\
	Grid_SetFlaggingField.o				\
	Grid_SetFlaggingFieldStaticRegions.o		\
	Grid_SetMinimumSupport.o			\
	Grid_ShockTubeInitializeGrid.o			\
	Grid_SolveForPotential.o			\
	Grid_SortParticlesByNumber.o			\
	Grid_SphericalInfallGetProfile.o		\
	Grid_SphericalInfallInitializeGrid.o		\
	Grid_StarParticleHandler.o			\
	Grid_SupernovaRestartInitialize.o		\
	Grid_TestGravityCheckResults.o			\
	Grid_TestGravitySphereCheckResults.o		\
	Grid_TestGravityInitializeGrid.o		\
	Grid_TestGravityMotionInitializeGrid.o		\
	Grid_TestGravitySphereInitializeGrid.o		\
	Grid_TestOrbitInitializeGrid.o			\
	Grid_TracerParticleCreateParticles.o		\
	Grid_TracerParticleOutputData.o			\
	Grid_TracerParticleSetVelocity.o		\
	Grid_WriteCube.o				\
	Grid_WriteCubeInterpolate.o			\
	Grid_WriteGridInterpolate.o			\
	Grid_Group_WriteGridInterpolate.o		\
	Grid_WriteMemoryMap.o				\
	Grid_WriteTaskMap.o				\
	Grid_WriteNewMovieData.o			\
	Grid_WriteStuff.o				\
	Grid_ZeldovichPancakeInitializeGrid.o		\
	Grid_ZeroSolutionUnderSubgrid.o			\
	Grid_constructor.o				\
	Grid_destructor.o				\
	Hierarchy_DeleteHierarchyEntry.o		\
	IdentifyNewSubgridsBySignature.o		\
	ImplosionInitialize.o				\
	InitializeCloudyCooling.o			\
	InitializeCosmicRayData.o			\
	InitializeEquilibriumCoolData.o			\
	InitializeGadgetEquilibriumCoolData.o		\
	InitializeLocal.o				\
	InitializeNew.o					\
	InitializeRadiationFieldData.o			\
	InitializeRateData.o				\
	InterpretCommandLine.o				\
	Isdigit.o					\
	KHInitialize.o					\
	LevelHierarchy_AddLevel.o			\
	ListIO.o					\
	MemoryAllocationRoutines.o			\
	MersenneTwister.o				\
	MultigridSolver.o				\
	NohInitialize.o					\
	OutputLevelInformation.o			\
	FastSib_PrepareGravitatingMassField.o		\
	PressurelessCollapseInitialize.o		\
	ProcMap.o					\
	ProtostellarCollapseInitialize.o		\
	ProtoSubgrid_AcceptableGrid.o			\
	ProtoSubgrid_CleanUp.o				\
	ProtoSubgrid_ComputeSecondDerivative.o		\
	ProtoSubgrid_ComputeSignature.o			\
	ProtoSubgrid_CopyFlaggedZonesFromGrid.o		\
	ProtoSubgrid_CopyToNewSubgrid.o			\
	ProtoSubgrid_FindGridsByZeroSignature.o		\
	ProtoSubgrid_ReturnNthLongestDimension.o	\
	ProtoSubgrid_ShrinkToMinimumSize.o		\
	ProtoSubgrid_constructor.o			\
	QuickSortAndDrag.o				\
	RadiatingShockInitialize.o			\
	RadiationFieldCalculateRates.o			\
	RadiationFieldUpdate.o				\
	ReadAllData.o					\
	Group_ReadAllData.o				\
	ReadDataHierarchy.o				\
	Group_ReadDataHierarchy.o			\
	ReadParameterFile.o				\
	ReadRadiationData.o				\
	ReadStarParticleData.o				\
	ReadUnits.o					\
	FastSib_RebuildHierarchy.o			\
	ReduceFragmentation.o				\
	ReportMemoryUsage.o				\
	RotatingCylinderInitialize.o			\
	SedovBlastInitialize.o				\
	SetDefaultGlobalValues.o			\
	SetDumpFlagsMask.o				\
	SetDumpFlagsOn.o				\
	ShockInABoxInitialize.o				\
	ShockPoolInitialize.o				\
	ShockTubeInitialize.o				\
	SphericalInfallInitialize.o			\
	SupernovaRestartInitialize.o			\
	SysMkdir.o					\
	TestGravityCheckResults.o			\
	TestGravityInitialize.o				\
	TestGravityMotion.o				\
	TestGravitySphereInitialize.o			\
	TestGravitySphereCheckResults.o			\
	TestOrbitInitialize.o				\
	TracerParticleCreation.o			\
	TurbulenceSimulationInitialize.o		\
	OptimizedUnigridCommunicationTranspose.o	\
	UpdateParticlePositions.o			\
	UpdateSymbolicLink.o				\
	WavePoolInitialize.o				\
	WriteAllData.o					\
	Group_WriteAllData.o				\
	WriteAllDataCubes.o				\
	WriteDataCubes.o				\
	WriteDataHierarchy.o				\
	Group_WriteDataHierarchy.o			\
	WriteHierarchyStuff.o				\
	WriteMemoryMap.o				\
	WriteMovieData.o				\
	WriteParameterFile.o				\
	WriteRadiationData.o				\
	WriteStarParticleData.o				\
	WriteStreamData.o				\
	WriteTaskMap.o					\
	WriteTracerParticleData.o			\
	WriteUnits.o					\
	ZeldovichPancakeInitialize.o			\
	E_ColumnFormat.o				\
	F_ColumnFormat.o				\
	P_ColumnFormat.o				\
	I_ColumnFormat.o				\
	arcsinh.o					\
	arccosh.o					\
	flow.o						\
	c_message.o					\
	f_message.o					\
	mused.o						\
	performance.o					\
	Options.o					\
	ENZO_Main.o					\
        calc_dt.o                                       \
	calc_photo_rates.o				\
	calc_rad.o					\
	calc_rates.o					\
        calcdiss.o					\
	chtable.o					\
	cicinterp.o					\
	cic_deposit.o					\
	cic_flag.o					\
	cic_interp.o					\
	colh2diss.o					\
	coll_rates.o					\
	comp_accel.o					\
	cool1d_cloudy.o					\
	cool_multi_time.o				\
	cool_time.o					\
        euler_nocheck.o					\
        expand_terms.o                                  \
	fortio.o					\
	fft66.o						\
	ffte4X.o					\
	prefort2.o					\
	cray_x1_3d_fft.o				\
	essl_3d_fft.o					\
	wrapper1d.o					\
	acml_1d_fft.o					\
	cray_x1_1d_fft.o				\
	essl_1d_fft.o					\
	ffte_1d_fft.o					\
	mkl_1d_fft.o					\
	s66_1d_fft.o					\
	s90_1d_fft.o					\
	sgi_1d_fft.o					\
	grid_cic.o					\
	int_lin3d.o					\
	int_spline.o					\
	interp1d.o					\
	interp2d.o					\
	interp3d.o					\
	interpolate.o					\
        inteuler.o					\
        intlgrg.o                                       \
        intrmp.o                                        \
        intvar.o                                        \
        lgrg.o                                          \
	mcooling.o					\
	mg_calc_defect.o				\
	mg_prolong.o					\
	mg_prolong2.o					\
	mg_relax.o					\
	mg_restrict.o					\
	mpi_error_file.o				\
	stop_all_cpus.o					\
	ngpinterp.o					\
        pgas2d.o					\
        pgas2d_dual.o					\
        ppm_lr.o					\
	project.o					\
	projplane.o					\
	prolong.o					\
	prolong_tsc.o					\
        remap.o                                         \
	smooth.o					\
	smooth_deposit.o				\
	star_maker1.o					\
	star_maker2.o					\
	star_maker3.o					\
	star_maker4.o					\
	star_maker5.o					\
	star_maker10.o					\
        tscint1d.o					\
        tscint2d.o					\
        tscint3d.o					\
        twoshock.o					\
        xeuler_sweep.o					\
        xlagrange_sweep.o                               \
        yeuler_sweep.o					\
        ylagrange_sweep.o                               \
        zeuler_sweep.o                                  \
	zeus_main.o					\
	zeus_utilities.o				\
	zeus_source.o					\
	zeus_xtransport.o				\
	zeus_ytransport.o				\
	zeus_ztransport.o				\
        zlagrange_sweep.o 				\
	FreeRealMem.o

H5OBJ =	SCB_ExternalBoundary_ReadExternalBoundary.o	\
	SCB_ExternalBoundary_WriteExternalBoundary.o	\
	write_bt.o	\
	write_bv.o	\
	read_bt.o	\
	read_bv.o	\
	Grid_CosmologySimulationInitializeGrid.o	\
	Grid_NestedCosmologySimulationInitializeGrid.o	\
	Grid_TurbulenceSimulationInitialize.o		\
	Grid_ReadGrid.o					\
	Grid_Group_ReadGrid.o				\
	Grid_WriteGrid.o				\
	Grid_Group_WriteGrid.o				\
	Grid_WriteGridX.o				\
	Grid_ReadRandomForcingFields.o			\
	ExtractSection.o				\
	OutputAsParticleData.o				\
	ProjectToPlane.o				\
	ReadAttr.o					\
	ReadFile.o					\
	ReadGridFile.o					\
	ReadIntFile.o					\
	WriteStringAttr.o

SROBJ = SRB_FileMover.o					\
	SRB_MakeDir.o					\
	SRB_Mkdir.o					\
	SRB_Clean.o					\
	SRB_Delete.o					\
	SRB_Mover.o					\
	SRB_Put.o

SDOBJ = SRB_DummyMakeDir.o				\
	SRB_DummyFileMover.o

RICKW = AnalysisBaseClass.o				\
	AnalysisBaseClass_HDF5Utils.o			\
	Grid_CreateFieldArray.o				\
	Grid_FlagGridArray.o				\
	Grid_FlagRefinedCells.o				\
	ealFloat.o					\
	ealInt.o

RICKT = analysisbaseclasstests.o			\
	ealfloattests.o					\
	ealinttests.o					\
	gridtests.o					\
	readgridtests.o

RAD =	\
	BlockSolve.o					\
	CalcEmiss.o					\
	CheckEmissivity.o				\
	ClearEmissivity.o				\
	CosmoIonizationInitialize.o			\
	DetermineParallelism.o				\
	DualCosmoIonizationInitialize.o 		\
	DualFLD_ComputeRadiationIntegrals.o		\
	DualFLD_ComputeOpacity.o			\
	DualFLD_ComputeTimeStep.o			\
	DualFLD_CrossSections.o 			\
	DualFLD_EnforceBoundary.o			\
	DualFLD_Evolve.o				\
	DualFLD_FillRates.o				\
	DualFLD_FInterface.o				\
	DualFLD_Initialize.o				\
	DualFLD_RadiationSource.o			\
	DualFLD_RadiationSpectrum.o			\
	DualFLD_SetupBoundary.o 			\
	DualFLD_SetupSystem.o				\
	DualFLD_WriteParameters.o			\
	DualFLD_constructor.o				\
	DualFLD_destructor.o                            \
	DualRadConstTestInitialize.o    		\
	DualRadStreamTestInitialize.o   		\
	DualRHIonizationTestInitialize.o		\
	EnzoVector.o					\
	EnzoVector_Exchange.o				\
	EnzoVector_ExchangeEnd.o			\
	EnzoVector_ExchangeStart.o			\
	EnzoVector_test.o				\
	FSMultiSourceInitialize.o			\
	FSProb_Dump.o					\
	FSProb_Evolve.o					\
	FSProb_FInterface.o				\
	FSProb_InitialGuess.o				\
	FSProb_Initialize.o				\
	FSProb_RadiationSource.o			\
	FSProb_SetupBoundary.o				\
	FSProb_SetupSystem.o				\
	FSProb_WriteParameters.o			\
	FSProb_constructor.o				\
	FSProb_destructor.o				\
	Grid_AccessBaryonFields.o			\
	Grid_CosmoIonizationInitializeGrid.o		\
	Grid_DualCosmoIonizationInitializeGrid.o	\
	Grid_DualRadConstTestInitializeGrid.o   	\
	Grid_DualRadStreamTestInitializeGrid.o   	\
	Grid_DualRHIonizationTestInitializeGrid.o	\
	Grid_FSMultiSourceInitializeGrid.o		\
	Grid_RHIonizationClumpInitializeGrid.o		\
	Grid_RHIonizationSteepInitializeGrid.o		\
	Grid_RHIonizationTestInitializeGrid.o		\
	Grid_RadHydroConstTestInitializeGrid.o		\
	Grid_RadHydroGreyMarshakWaveInitializeGrid.o	\
	Grid_RadHydroPulseTestInitializeGrid.o		\
	Grid_RadHydroRadShockInitializeGrid.o		\
	Grid_RadHydroStreamTestInitializeGrid.o		\
	InexactNewton.o					\
	InexactNewton_InexactNewtonForce.o		\
	InexactNewton_LinesearchStepSize.o		\
	InexactNewton_Solve.o				\
	InexactNewton_constructor.o			\
	InexactNewton_destructor.o			\
	NullProblem.o					\
	RHIonizationClumpInitialize.o			\
	RHIonizationSteepInitialize.o			\
	RHIonizationTestInitialize.o			\
	RadHydroConstTestInitialize.o			\
	RadHydroGreyMarshakWaveInitialize.o		\
	RadHydroPulseTestInitialize.o			\
	RadHydroRadShockInitialize.o			\
	RadHydroStreamTestInitialize.o			\
	RadiationGetUnits.o				\
	UnigridClearEmissivity.o			\
	gFLDProblem_AnalyticalEqns.o			\
	gFLDProblem_ChemistrySource.o			\
	gFLDProblem_ComputeRHS.o			\
	gFLDProblem_ComputeRadiationIntegrals.o		\
	gFLDProblem_ComputeTemperature.o		\
	gFLDProblem_ComputeTimeStep.o			\
	gFLDProblem_CrossSections.o			\
	gFLDProblem_DiffRHS.o				\
	gFLDProblem_Dump.o				\
	gFLDProblem_EnforceBoundary.o			\
	gFLDProblem_Evolve.o				\
	gFLDProblem_FInterface.o			\
	gFLDProblem_GasEnergySource.o			\
	gFLDProblem_InitialGuess.o			\
	gFLDProblem_Initialize.o			\
	gFLDProblem_LocRHS.o				\
	gFLDProblem_LocalJac.o				\
	gFLDProblem_LocalRHS.o				\
	gFLDProblem_MatrixEntries.o			\
	gFLDProblem_Opacity.o				\
	gFLDProblem_RadiationSource.o			\
	gFLDProblem_RadiationSpectrum.o			\
	gFLDProblem_SetNewtonBCs.o			\
	gFLDProblem_SetupBoundary.o			\
	gFLDProblem_UpdateBoundary.o			\
	gFLDProblem_WriteParameters.o			\
	gFLDProblem_constructor.o			\
	gFLDProblem_destructor.o			\
	gFLDProblem_lsetup.o				\
	gFLDProblem_lsolve.o				\
	gFLDProblem_nlresid.o				\
	gFLDSplit_AnalyticChemistry.o			\
	gFLDSplit_ChemistrySource.o			\
	gFLDSplit_ComputeRadiationIntegrals.o		\
	gFLDSplit_ComputeTemperature.o			\
	gFLDSplit_ComputeTimeStep.o			\
	gFLDSplit_CrossSections.o			\
	gFLDSplit_Dump.o				\
	gFLDSplit_EnforceBoundary.o			\
	gFLDSplit_Evolve.o				\
	gFLDSplit_FillRates.o				\
	gFLDSplit_FInterface.o				\
	gFLDSplit_GasEnergySource.o			\
	gFLDSplit_InitialGuess.o			\
	gFLDSplit_Initialize.o				\
	gFLDSplit_Opacity.o				\
	gFLDSplit_RadiationSource.o			\
	gFLDSplit_RadiationSpectrum.o			\
	gFLDSplit_SetupBoundary.o			\
	gFLDSplit_SetupSystem.o				\
	gFLDSplit_WriteParameters.o			\
	gFLDSplit_constructor.o				\
	gFLDSplit_destructor.o

NEW_THREAD = \
	NEW_cool1d.o					\
	NEW_cool1d_multi.o

SPARE = \
	NEW_multi_cool.o


OMP_CRIT = \
	OMP_CRITICAL_Grid_ComputeCoolingTime.o		\
	OMP_CRITICAL_Grid_SolveRadiativeCooling.o	\
	OMP_CRITICAL_Grid_SolveRateAndCoolEquations.o	\
	OMP_CRITICAL_Grid_SolveRateEquations.o

OMP_BASE = \
	OMP_NEW_multi_cool.o				\
	OMP_solve_rate_cool.o				\
	OMP_FastSib_EvolveLevelRoutinesOptimized.o	\
	OMP_Grid_CopyBaryonFieldToOldBaryonField.o	\
	OMP_Grid_GetProjectedBoundaryFluxes.o		\
	OMP_Grid_ProjectSolutionToParentGrid.o		\
	OMP_Grid_UpdateParticlePosition.o		\
	OMP_Grid_UpdateParticleVelocity.o		\
	OMP_fft_utils.o					\
	OMP_rotate2d.o					\
	OMP_rotate3d.o					\
	OMP_solve_cool.o				\
	OMP_solve_rate.o				\
	OMP_utilities.o					\
	OMP_wrapper2d.o					\
	OMP_wrapper3d.o

OMP_X = \
	OMP_Grid_CorrectForRefinedFluxes.o		\
	OMP_FluxFix_Grid_CorrectForRefinedFluxes.o	\
	OMP_EvolveLevel.o

OMP_FORT = \
	FORTRAN_Grid_SolveHydroEquations.o		\
	OMP_local1_ppm_de.o

OMP_WISE2 = \
	WISE_Grid_SolveHydroEquations.o			\
	OMP_Grid_SolvePPM_DE.o				\
	jw_euler.o					\
	jw_intvar.o					\
	jw_inteuler.o					\
	flux_hllc.o					\
	flux_hll.o					\
	calc_eigen.o					\
	flux_twoshock.o					\
	intprim.o					\
	intpos.o					\
	Grid_xEulerSweep_W2.o				\
	Grid_yEulerSweep_W2.o				\
	Grid_zEulerSweep_W2.o

OMP_WISE1 = \
	WISE_Grid_SolveHydroEquations.o			\
	OMP_Grid_SolvePPM_DE.o				\
	Grid_xEulerSweep_W1.o				\
	Grid_yEulerSweep_W1.o				\
	Grid_zEulerSweep_W1.o

OLD_THREAD = \
	cool1d.o					\
	cool1d_multi.o					\
	multi_cool.o

OLD_CRIT = \
	Grid_ComputeCoolingTime.o			\
	Grid_SolveRadiativeCooling.o			\
	Grid_SolveRateAndCoolEquations.o		\
	Grid_SolveRateEquations.o

OLD_BASE = \
	FastSib_EvolveLevelRoutinesOptimized.o		\
	Grid_CopyBaryonFieldToOldBaryonField.o		\
	Grid_GetProjectedBoundaryFluxes.o		\
	Grid_ProjectSolutionToParentGrid.o		\
	Grid_UpdateParticlePosition.o			\
	Grid_UpdateParticleVelocity.o			\
	fft_utils.o					\
	rotate2d.o					\
	rotate3d.o					\
	solve_rate_cool.o				\
	solve_cool.o					\
	solve_rate.o					\
	utilities.o					\
	wrapper2d.o					\
	wrapper3d.o

OLD_X = \
	Grid_CorrectForRefinedFluxes.o			\
	FluxFix_Grid_CorrectForRefinedFluxes.o		\
	EvolveLevel.o

OLD_WISE2 = \
	WISE_Grid_SolveHydroEquations.o			\
	Grid_SolvePPM_DE.o				\
	jw_euler.o					\
	jw_intvar.o					\
	jw_inteuler.o					\
	flux_hllc.o					\
	flux_hll.o					\
	calc_eigen.o					\
	flux_twoshock.o					\
	intprim.o					\
	intpos.o					\
	Grid_xEulerSweep_W2.o				\
	Grid_yEulerSweep_W2.o				\
	Grid_zEulerSweep_W2.o

OLD_WISE1 = \
	WISE_Grid_SolveHydroEquations.o			\
	Grid_SolvePPM_DE.o				\
	Grid_xEulerSweep_W1.o				\
	Grid_yEulerSweep_W1.o				\
	Grid_zEulerSweep_W1.o

OLD_FORT = \
	FORTRAN_Grid_SolveHydroEquations.o		\
	ppm_de.o

OLD_WOC = \
	WOC_Grid_SolveHydroEquations.o			\
	Grid_PPMDirectEuler.o				\
	feuler_sweep.o

ifeq ($(HYDRO),FORTRAN)
OLD_HYDRO := ${OLD_FORT}
OMP_HYDRO := ${OMP_FORT}
endif

ifeq ($(HYDRO),WISE1)
OLD_HYDRO := ${OLD_WISE1}
OMP_HYDRO := ${OMP_WISE1}
endif

ifeq ($(HYDRO),WISE2)
OLD_HYDRO := ${OLD_WISE2}
OMP_HYDRO := ${OMP_WISE2}
endif

ifeq ($(ENZO_INT),SMALL_INTS)
  ifeq ($(ENZO_AMR),yes)
    EXENAME := enzo5a_32
  else
    EXENAME := enzo5u_32
  endif
else
  ifeq ($(ENZO_AMR),yes)
    EXENAME := enzo5a_64
  else
    EXENAME := enzo5u_64
  endif
endif

ifeq ($(ENZO_RHD),yes)
ifeq ($(ENZO_MODE),hybrid)
  OBJ     := ${F90X} ${XOBJ} ${H5OBJ} ${SDOBJ} ${RICKW} ${RAD} ${NEW_THREAD} ${OMP_CRIT} ${OMP_BASE} ${OMP_X} ${OMP_HYDRO}
else
  OBJ     := ${F90X} ${XOBJ} ${H5OBJ} ${SDOBJ} ${RICKW} ${RAD} ${OLD_THREAD} ${OLD_CRIT} ${OLD_BASE} ${OLD_X} ${OLD_HYDRO}
endif
else
ifeq ($(ENZO_MODE),hybrid)
  OBJ     := ${F90X} ${XOBJ} ${H5OBJ} ${SDOBJ} ${RICKW} ${NEW_THREAD} ${OMP_CRIT} ${OMP_BASE} ${OMP_X} ${OMP_HYDRO}
else
  OBJ     := ${F90X} ${XOBJ} ${H5OBJ} ${SDOBJ} ${RICKW} ${OLD_THREAD} ${OLD_CRIT} ${OLD_BASE} ${OLD_X} ${OLD_HYDRO}
endif
endif


#-------------------  macro definitions --------------------------------

.SUFFIXES:
.SUFFIXES: .src .src90 .c .f .o .C

#HEADER  = CosmologyParameters.h ExternalBoundary.h Fluxes.h \
#	  Grid.h GridList.h Hierarchy.h LevelHierarchy.h ProtoSubgrid.h \
#	  TopGridData.h global_data.h macros_and_parameters.h typedefs.h region.h \
#	  CoolData.h RateData.h StarParticleData.h RadiationFieldData.h \
#	  fortran.def

HEADER  = \
	CoolData.h \
	CosmologyParameters.h \
	ExternalBoundary.h \
	Fluxes.h \
	fortran.def \
	global_data.h \
	Grid.h \
	GridList.h \
	Hierarchy.h \
	LevelHierarchy.h \
	macros_and_parameters.h\
	ProtoSubgrid.h \
	RadiationFieldData.h \
	RateData.h\
	region.h \
	StarParticleData.h \
	TopGridData.h \
	typedefs.h

#------------------ include machine-dependent macros -------------------

SYS_INCLUDE = ${ENZO_SYSTEM}_${ENZO_COMPILER}.${ENZO_MODE}${ENZO_SPEC}
MAKEFILE_INCLUDE = Macros/Make.${SYS_INCLUDE}
include ${MAKEFILE_INCLUDE}

#--------------------  implicit rules ----------------------------------

#include ${MAKEFILE_RULES}

ENZO_PARAM := -D__max_subgrids=${MAX_GRIDS} -D__max_baryons=${MAX_BARYONS} -D__max_colors=${MAX_COLORS} -D__max_cube_plus_7=${MAX_CUBE_PLUS_7} -D__max_cpu_per_node=${MAX_CPU_PER_NODE}

.src90.o:
	${CPP} ${RHD_OPT} ${FC_DEF} -D${FFT_OPT} -D${ENZO_PPR} -D${ENZO_XPR} -P ${CPP_FLAG} ${ENZO_PARAM} $< > $*.f90
	${FC} ${F90SRC} ${FC_OPT} ${FREE90} ${DEBUG} ${OPTS} $*.f90
.src.o:
	${CPP} ${RHD_OPT} ${FC_DEF} -D${FFT_OPT} -D${ENZO_PPR} -D${ENZO_XPR} -P ${CPP_FLAG} ${ENZO_PARAM} $< > $*.f
	${FC} ${FC_OPT} ${DEBUG} ${OPTS} $*.f 
.C.o:
	${CCxx} ${CCxx_DEF} ${CCxx_OPT} ${AMR_OPT} ${RHD_OPT} -D${OOC_OPT} -D${UXP_OPT} -D${ENZO_INT} -D${INIT_PPR} -D${ENZO_PPR} -D${ENZO_XPR} -D${ENZO_END} -D${ENZO_TAG} -D${TASK_MAP} -D${IO_PR} -D${IO_PACK} -D${IO_OPT} -D${IO_MEM_BUFF_OUT} -D${IO_MEM_BUFF_INP} -D${ISO_GRAV} -D${JB_PERF} ${ENZO_PARAM} ${DEBUG} ${OPTS} $*.C
.c.o:
	${CC} ${CC_DEF} ${CC_OPT} -D${ENZO_INT} -D${ENZO_PPR} -D${ENZO_XPR} -D${ENZO_END} -D${ENZO_TAG} -D${TASK_MAP} -D${USE_LOADL} ${DEBUG} ${OPTS} -c $<

#-------------------------  main target  -------------------------------

${EXENAME}: ${OBJ}
	${LDR} ${LDR_OPT} ${DEBUG} -o ${EXENAME} ${OBJ} ${LIB}
	ln -fs ${EXENAME} enzo.exe


#-------------------------  documentation ------------------------------

Doc:
	@echo "VERSION   : " ${VERSION}          > Compile_Options
	@echo "DATE      : " ${DATE}            >> Compile_Options
	@echo "CWD       : " ${CWD}             >> Compile_Options
	@echo "SYSTEM    : " ${ENZO_SYSTEM}     >> Compile_Options
	@echo "COMPILER  : " ${ENZO_COMPILER}   >> Compile_Options
	@echo " "                               >> Compile_Options
	@echo "ENZO_MODE : " ${ENZO_MODE}       >> Compile_Options
	@echo "ENZO_SPEC : " ${ENZO_SPEC}       >> Compile_Options
	@echo "ENZO_INTS : " ${ENZO_INT}        >> Compile_Options
	@echo "ENZO_TAGS : " ${ENZO_TAG}        >> Compile_Options
	@echo "INIT_PPR  : " ${INIT_PPR}        >> Compile_Options
	@echo "ENZO_PPR  : " ${ENZO_PPR}        >> Compile_Options
	@echo "ENZO_XPR  : " ${ENZO_XPR}        >> Compile_Options
	@echo "ENZO_AMR  : " ${ENZO_AMR}        >> Compile_Options
	@echo "ENZO_RHD  : " ${ENZO_RHD}        >> Compile_Options
	@echo "PACK_AMR  : " ${PACK_AMR}        >> Compile_Options
	@echo "PREC_IO   : " ${PREC_IO}         >> Compile_Options
	@echo "TASK_MAP  : " ${IBM_TASKMAP}     >> Compile_Options
	@echo "ENZO_SRB  : " ${ENZO_SRB}        >> Compile_Options
	@echo "ENZO_END  : " ${ENZO_END}        >> Compile_Options
	@echo "ENZO_PAR  : " ${ENZO_PARAM}      >> Compile_Options
	@echo " "                               >> Compile_Options
	@echo "RHD_OPT   : " ${RHD_OPT}         >> Compile_Options
	@echo "AMR_OPT   : " ${AMR_OPT}         >> Compile_Options
	@echo "OOC_OPT   : " ${OOC_OPT}         >> Compile_Options
	@echo "UXP_OPT   : " ${UXP_OPT}         >> Compile_Options
	@echo "FFT_OPT   : " ${FFT_OPT}         >> Compile_Options
	@echo " "                               >> Compile_Options
	@echo "INCLUDE   : " ${MAKEFILE_INCLUDE}  >> Compile_Options
	@echo " "                               >> Compile_Options
	@echo "CPP       : " ${CPP}             >> Compile_Options
	@echo "CPP_FLAG  : " ${CPP_FLAG}        >> Compile_Options
	@echo "OPTIONS   : " ${OPTS}            >> Compile_Options
	@echo " "                               >> Compile_Options
	@echo "FC        : " ${FC}              >> Compile_Options
	@echo "FC_OPT    : " ${FC_OPT}          >> Compile_Options
	@echo "FC_DEF    : " ${FC_DEF}          >> Compile_Options
	@echo "FREE90    : " ${FREE90}          >> Compile_Options
	@echo " "                               >> Compile_Options
	@echo "CC        : " ${CC}              >> Compile_Options
	@echo "CC_OPT    : " ${CC_OPT}          >> Compile_Options
	@echo "CC_DEF    : " ${CC_DEF}          >> Compile_Options
	@echo " "                               >> Compile_Options
	@echo "C++       : " ${CCxx}            >> Compile_Options
	@echo "C++_OPT   : " ${CCxx_OPT}        >> Compile_Options
	@echo "C++_DEF   : " ${CCxx_DEF}        >> Compile_Options
	@echo " "                               >> Compile_Options
	@echo "LDR       : " ${LDR}             >> Compile_Options
	@echo "LDR_OPT   : " ${LDR_OPT}         >> Compile_Options
	@echo "LIB       : " ${LIB}             >> Compile_Options
	@echo " "                               >> Compile_Options
	@echo "HDF5      : " ${HDF5_INSTALL}    >> Compile_Options
	@echo "HDF5_PP   : " ${HDF5_INSTALL_PP} >> Compile_Options
	@echo "HYPRE     : " ${HYPRE_INSTALL}   >> Compile_Options

#---------------------  special compilation rules ----------------------


#----------------------------- help ------------------------------------

help:
	@echo Type 'make' to generate executable
	@echo Type 'make clean' to remove object files, etc

#----------------------------- clean -----------------------------------

clean:
	rm -f *.o *.mod *.f *.f90 *~ ${EXENAME} Compile_Options load_map

#-----------------------------------------------------------------------

${OBJ}: ${HEADER}
