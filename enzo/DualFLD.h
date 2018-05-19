/*****************************************************************************
 * Copyright 2013 Daniel R. Reynolds                                         *
 *                                                                           *
 * This software is released under the terms of the "Enzo Public License"    *
 * in the accompanying LICENSE file.                                         *
 *****************************************************************************/
/***********************************************************************
   Two-Group (X-ray + UV) Flux-Limited Diffusion Implicit Problem Class

   author: Daniel R. Reynolds
   date:   January 2013

   PURPOSE: This class defines problem-specific data and functions 
            for an implicit, 2-group flux-limited diffusion solves.

            The variables are stored in the following order: 
               0 -> UV radiation energy density
               1 -> X-ray radiation energy density
************************************************************************/
#ifdef RAD_HYDRO
#ifndef DUAL_FLD_PROBLEM_DEFINED__
#define DUAL_FLD_PROBLEM_DEFINED__

#include "ImplicitProblemABC_preincludes.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "Hierarchy.h"
#include "TopGridData.h"
#include "EnzoVector.h"
#include "ImplicitProblemABC.h"


class DualFLD : public virtual ImplicitProblemABC {

 private:
  
  // overall time spent in solver and components
  float RTtime;
  float HYPREtime;
  
  // HYPRE Struct-specific data
  Eint32 mattype;                // HYPRE matrix type for solve
  Eint32 stSize;                 // stencil size
  HYPRE_StructGrid grid;         // HYPRE grid object for setup
  HYPRE_StructStencil stencil;   // stencil object

  // HYPRE Solver-specific data
  float  sol_tolerance_Xr;       // desired solver tolerance
  float  sol_tolerance_UV;
  Eint32 sol_MGmaxit_Xr;         // maximum number of MG iterations
  Eint32 sol_MGmaxit_UV;
  Eint32 sol_PCGmaxit_Xr;        // maximum number of PCG iterations
  Eint32 sol_PCGmaxit_UV;
  Eint32 sol_rlxtype_Xr;         // relaxation type:
                                 //    0,1 -> weighted Jacobi
                                 //    2,3 -> red-black Gauss-Seidel
  Eint32 sol_rlxtype_UV;
  Eint32 sol_npre_Xr;            // num. pre-relaxation sweeps
  Eint32 sol_npre_UV;
  Eint32 sol_npost_Xr;           // num. post-relaxation sweeps
  Eint32 sol_npost_UV;
  Eint32 sol_printl;             // print output level
  Eint32 sol_log;                // amount of logging
  Eint32 SolvIndices[3][2];      // L/R edge indices of subdomain in global mesh
                                 // Note: these INCLUDE Dirichlet zones, even 
                                 //   though those are not included as active 
                                 //   data in the vectors or physics routines.
  int SolvOff[3];                // offset between HYPRE mesh and active mesh; 
                                 //   typically 0, but will be 1 for inclusion 
                                 //   of Dirichlet zones in HYPRE grid.

  // HYPRE interface temporary data
  HYPRE_StructMatrix P;          // holds radiation matrix
  HYPRE_StructVector rhsvec;     // holds radiation rhs vector
  HYPRE_StructVector solvec;     // holds radiation solution vector
  Eflt64 *matentries;            // holds radiation matrix entries
  Eflt64 *rhsentries;            // linear system rhs entries
  Eflt64 *HYPREbuff;             // holds contiguous sections of rhs/sol

  // HYPRE solver diagnostics
  int totIters;                  // total MG iterations for solves

  // General problem grid information
  bool OnBdry[3][2]; // denotes if proc owns piece of boundary
  int rank;          // Rank of self-gravity problem
  int layout[3];     // number of procs in each dim (1-based)
  int location[3];   // location of this proc in each dim (0-based)
  int NBors[3][2];   // process IDs of L/R neighbors in each dim
  int LocDims[3];    // implicit problem local dims (no ghost or bdry cells)
  int ArrDims[3];    // local array sizes (includes ghost and bdry cells)
  int GhDims[3][2];  // ghost cells at each face (includes Dirichlet bdry zones)
  int GlobDims[3];   // implicit problem global dimensions (active cells only)
  float dx[3];             // mesh size in each dimension
  float EdgeVals[3][2];    // L/R edges of this proc's subdomain
  float *UVBdryVals[3][2];   // boundary values for UV BCs
  float *XrBdryVals[3][2];   // boundary values for Xray BCs

  // time-stepping related data
  float initdt;          // initial radiation time step size
  float maxdt;           // maximum radiation/chemistry/heating time step size
  float mindt;           // minimum radiation/chemistry/heating time step size
  float dtfac[2];        // desired relative change in fields per step
  float dtnorm;          // norm choice for computing relative change:
                         //    0 -> max pointwise norm (default)
                         //   >0 -> rms p-norm over entire domain
  float dtgrowth;        // time step growth factor (1 < dtgrowth < 10)
  float tnew;            // new time
  float told;            // old time
  float dt;              // time step size
  float dtrad;           // radiation time step size (subcycled)
  float theta;           // implicitness parameter (1->BE, 0.5->CN, 0->FE)
  EnzoVector *sol;       // solution vector
  EnzoVector *U0;        // old time-level state
  EnzoVector *extsrc;    // temporary vector holding external forcing sources
  EnzoVector *opacity;
  
  // problem defining data
  bool XrayOnly;       // flag to enable/disable UV radiation
  bool XrStatic;       // flag to denote a static Xray radiation field
  bool UVStatic;       // flag to denote a static UV radiation field
  int Nchem;           // number of chemical species (non-negative integer)
  float NGammaDotUV;   // UV source strength (photons/sec)
  float NGammaDotXr;   // X-ray source strength (photons/sec)
  float EtaRadius;     // ionization source radius
  float EtaCenter[3];  // ionization source location

  // cosmology and scaling constants
  FLOAT a;             // cosmology expansion coefficient
  FLOAT a0;            // cosmology expansion coefficient (old time)
  FLOAT adot;          // time-derivative of a
  FLOAT adot0;         // time-derivative of a (old time)
  float aUnits;        // expansion parameter scaling
  bool  autoScale;     // flag to enable/disable automatic scaling factors
  float UVScale;       // scaling factor for UV radiation energy density
  float UVUnits;       // UV radiation energy density unit conversion factor
  float UVUnits0;      // UV radiation energy density unit conversion factor
  float XrScale;       // scaling factor for X-ray radiation energy density
  float XrUnits;       // X-ray radiation energy density unit conversion factor
  float XrUnits0;      // X-ray radiation energy density unit conversion factor
  float NiUnits;       // species density unit conversion factor
  float NiUnits0;      // species density unit conversion factor

  float DenUnits;      // density scaling factor
  float LenUnits;      // length scaling factor
  float TimeUnits;     // time scaling factor
  float VelUnits;      // velocity scaling factor
  float DenUnits0;     // density scaling factor
  float LenUnits0;     // length scaling factor

  // chemistry constants
  float HFrac;         // Fraction of matter composed of Hydrogen

  // storage for integrals over radiation spectrum (set during initialization)
  float hnu0_HI;       // HI ionization threshold (eV)
  float hnu0_HeI;      // HeI ionization threshold (eV)
  float hnu0_HeII;     // HeII ionization threshold (eV)
  int UVSpectrum;      // UV radiation spectrum choice:
                       //   2 -> PopII SED (like Ricotti etal, ApJ 575:33–48, 2002)
                       //   1 -> 1e5 black body spectrum
                       //   0 -> simple power law spectrum (exponent -1.5)
                       //  -1 -> monochromatic spectrum at UVFrequency
  float UVFrequency;   // UV radiation frequency (if UVSpectrum = -1)
  int XrSpectrum;      // X-ray spectrum choice (same values as above)
  float XrFrequency;   // X-ray radiation frequency (if XrSpectrum = -1)
  float RadIntUV[7];   // UV spectrum integrals:
                       //   0: int_{nu0}^{inf} chiUV(nu) d nu
                       //   1: int_{nu0}^{inf} chiUV(nu)*sigmaHI(nu) dnu
                       //   2: int_{nu0}^{inf} chiUV(nu)*sigmaHI(nu)/nu dnu
                       //   3: int_{nu0}^{inf} chiUV(nu)*sigmaHeI(nu) dnu
                       //   4: int_{nu0}^{inf} chiUV(nu)*sigmaHeI(nu)/nu dnu
                       //   5: int_{nu0}^{inf} chiUV(nu)*sigmaHeII(nu) dnu
                       //   6: int_{nu0}^{inf} chiUV(nu)*sigmaHeII(nu)/nu dnu
  float RadIntXr[7];   // X-ray spectrum integrals (same form as above)

  // stored arrays for increased efficiency

  // private computation routines
  int   EnforceBoundary(EnzoVector *u);
  int   SetupSystem(float *rhsnorm, int UV_Xr);
  int   ComputeOpacity(float *HI, float *HeI, float *HeII);
  int   RadiationSource();
  float RadiationSpectrum(int ispectrum, float nu);
  float CrossSections(float nu, int species);
  int   ComputeRadiationIntegrals();
  int   FillRates(EnzoVector *u, float *rho, float *HI, float *HeI, 
		  float *HeII, float *phHI, float *phHeI, float *phHeII, 
		  float *PhotoGamma, float *dissH2I);
  int   RadStep(HierarchyEntry *ThisGrid, int eta_set, int UV_Xr);
  

 public:

  // boundary type in each dimension, face:
  //    0->periodic
  //    1->dirichlet
  //    2->neumann
  int XrBdryType[3][2];
  int UVBdryType[3][2];

  ///////////////////////////////////////
  // FLD-Specific Routines

  // Constructor
  DualFLD();
  
  // Destructor
  ~DualFLD();

  // Problem Initializer
  int Initialize(HierarchyEntry &TopGrid, TopGridData &MetaData);
  
  // Problem Evolver
  int Evolve(HierarchyEntry *ThisGrid, float deltat);
  
  // Write module parameters to file
  int WriteParameters(FILE *fptr);

  // Problem Boundary Condition setup (called once or at each time step, 
  //    must be called for each locally-owned external face separately)
  int SetupBoundary(int Dimension, int Face, int XrUV, int BdryConst, float *BdryData);

  // Return the maximum rad-hydro time step size
  float ComputeTimeStep(EnzoVector *uold, EnzoVector *unew, int flag);

};


#endif
#endif
