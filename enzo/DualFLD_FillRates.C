/*****************************************************************************
 * Copyright 2013 Daniel R. Reynolds                                         *
 *                                                                           *
 * This software is released under the terms of the "Enzo Public License"    *
 * in the accompanying LICENSE file.                                         *
 *****************************************************************************/
/***********************************************************************
   Two-Group (X-ray + UV) Flux-Limited Diffusion Implicit Problem Class
   FillRates routine.

   author: Daniel R. Reynolds
   date:   January 2013

   PURPOSE: Fills the photo-ionization and photo-heating arrays used by
            chemistry and cooling routines using time-averaged internal
            values for the radiation.
 
   NOTE: In order to save on memory, the photo-heating rates are 
         combined into a single rate, and scaled by the current number 
         density of HI, to be later unpacked by rescaling back with HI.  
         This loses accuracy in the case that during chemistry 
         subcycling the chemistry changes significantly, since we retain 
         the initial rate scaling but use updated HI values in rescaling.
************************************************************************/
#ifdef RAD_HYDRO
#include "DualFLD.h"


int DualFLD::FillRates(EnzoVector *u, float *rho, float *HI, float *HeI, 
		       float *HeII, float *phHI, float *phHeI, 
		       float *phHeII, float *photogamma, float *dissH2I) {

  // get local mesh description
  int usz[4], ghXl, ghXr, ghYl, ghYr, ghZl, ghZr;
  u->size(&usz[0], &usz[1], &usz[2], &usz[3], 
	  &ghXl, &ghXr, &ghYl, &ghYr, &ghZl, &ghZr);
  if (usz[0] != LocDims[0]) {
    fprintf(stderr,"FillRates error: x0 vector dims do not match\n");
    return FAIL;  }
  if (usz[1] != LocDims[1]) {
    fprintf(stderr,"FillRates error: x1 vector dims do not match\n");
    return FAIL;  }
  if (usz[2] != LocDims[2]) {
    fprintf(stderr,"FillRates error: x2 vector dims do not match\n");
    return FAIL;  }
  if ((usz[0]+ghXl+ghXr) != ArrDims[0]) {
    fprintf(stderr,"FillRates error: x0 vector sizes do not match\n");
    return FAIL;  }
  if ((usz[1]+ghYl+ghYr) != ArrDims[1]) {
    fprintf(stderr,"FillRates error: x1 vector sizes do not match\n");
    return FAIL;  }
  if ((usz[2]+ghZl+ghZr) != ArrDims[2]) {
    fprintf(stderr,"FillRates error: x2 vector sizes do not match\n");
    return FAIL;  }

  // set some physical constants
  float c        = 2.99792458e10;    // speed of light [cm/s]
  float hp       = 6.6260693e-27;    // Planck's constant [ergs*s]
  float mp       = 1.67262171e-24;   // mass of a proton [g]
  float ev2erg   = 1.60217653e-12;   // conversion constant from eV to ergs
  float dom      = DenUnits*a*a*a/mp;
  float tbase1   = TimeUnits;
  float xbase1   = LenUnits/a/aUnits;
  float dbase1   = DenUnits*a*a*a*aUnits*aUnits*aUnits;
  float coolunit = aUnits*aUnits*aUnits*aUnits*aUnits * xbase1*xbase1
                 * mp*mp / tbase1/tbase1/tbase1 / dbase1;
  float rtunits  = ev2erg/TimeUnits/coolunit/dom;
  float UVUn     = (UVUnits+UVUnits0)*0.5;
  float XrUn     = (XrUnits+XrUnits0)*0.5;


  // access radiation energy density array
  float *UV=NULL, *Xr=NULL;
  Xr = u->GetData(0);
  if (!XrayOnly)  UV = u->GetData(1);


  // photogamma array temporarily holds electron fraction for secondary ionization rates
  float *xe = photogamma;

  // compute the size of the fields
  int i, size=1;
  for (i=0; i<rank; i++)  size *= ArrDims[i];

  // fill electron fraction over grid
  if (Nchem == 1)
    for (i=0; i<size; i++)
      xe[i] = max(1.0 - HI[i]/(rho[i]*HFrac), 1e-4);   // xe = nHII / nH
  if (Nchem == 3)
    for (i=0; i<size; i++)   // xe = (nHII + nHeII/4 + nHeIII/2) / (nH + 2*nHe)
      xe[i] = max((rho[i]*(HFrac+1.0)*0.5 - HI[i] - 0.5*HeI[i] - 0.25*HeII[i])
		  / (rho[i]*(HFrac+1.0)*0.5), 1e-4);

  // fill HI photo-ionization rate
  float pHIconstUV = c*TimeUnits*UVUn*RadIntUV[2]/hp/RadIntUV[0];
  float pHIconstXr = c*TimeUnits*XrUn*RadIntXr[1]/RadIntXr[0]/(hnu0_HI*ev2erg);
  for (i=0; i<size; i++)
    phHI[i] = Xr[i]*pHIconstXr*0.3908*pow((1.0-pow(xe[i], 0.4092)), 1.7592);
  if (!XrayOnly)
    for (i=0; i<size; i++)  phHI[i] += UV[i]*pHIconstUV;

  // fill HeI and HeII photo-ionization rates
  float pHeIconstUV  = c*TimeUnits*UVUn*RadIntUV[4]/hp/RadIntUV[0];
  float pHeIIconstUV = c*TimeUnits*UVUn*RadIntUV[6]/hp/RadIntUV[0];
  float pHeIconstXr  = c*TimeUnits*XrUn*RadIntXr[3]/RadIntXr[0]/(hnu0_HeI*ev2erg);
  if (RadiativeTransferHydrogenOnly == FALSE) {
    for (i=0; i<size; i++)  
      phHeI[i] = Xr[i]*pHeIconstXr*0.0554*pow((1.0-pow(xe[i], 0.4614)), 1.6660);
    if (!XrayOnly) {
      for (i=0; i<size; i++)  phHeI[i] += UV[i]*pHeIconstUV;
      for (i=0; i<size; i++)  phHeII[i] = UV[i]*pHeIIconstUV;
    } else {
      for (i=0; i<size; i++)  phHeII[i] = 0.0;
    }
  }
   
  // fill photo-heating rate
  float phScaleUV    = c*TimeUnits*UVUn/RadIntUV[0]/VelUnits/VelUnits/mp/rtunits;
  float GHIconstUV   = phScaleUV*(RadIntUV[1] - 13.6*ev2erg/hp*RadIntUV[2]);
  float GHeIconstUV  = phScaleUV*(RadIntUV[3] - 24.6*ev2erg/hp*RadIntUV[4]);
  float GHeIIconstUV = phScaleUV*(RadIntUV[5] - 54.4*ev2erg/hp*RadIntUV[6]);
  float phScaleXr    = c*TimeUnits*XrUn/RadIntXr[0]/VelUnits/VelUnits/mp/rtunits*0.9971;
  float GHIconstXr   = phScaleXr*(RadIntXr[1] - 13.6*ev2erg/hp*RadIntXr[2]);
  float GHeIconstXr  = phScaleXr*(RadIntXr[3] - 24.6*ev2erg/hp*RadIntXr[4]);
  float GHeIIconstXr = phScaleXr*(RadIntXr[5] - 54.4*ev2erg/hp*RadIntXr[6]);
  if (Nchem == 1) {
    for (i=0; i<size; i++)
      photogamma[i] = Xr[i]*GHIconstXr*(1.0 - pow(1.0 - pow(xe[i], 0.2663), 1.3163));
    if (!XrayOnly) 
      for (i=0; i<size; i++)  photogamma[i] += UV[i]*GHIconstUV;
  }
  if (Nchem == 3) {
    for (i=0; i<size; i++)  photogamma[i] = 
	Xr[i]/HI[i]*(GHIconstXr*HI[i] + GHeIconstXr*HeI[i] + GHeIIconstXr*HeII[i])
		   *(1.0 - pow(1.0 - pow(xe[i], 0.2663), 1.3163));
    if (!XrayOnly)
      for (i=0; i<size; i++)  photogamma[i] += 
	UV[i]/HI[i]*(GHIconstUV*HI[i] + GHeIconstUV*HeI[i] + GHeIIconstUV*HeII[i]);
  }

  // fill H2 dissociation rate (zero for grey FLD problems)
  if (MultiSpecies > 1) 
    for (i=0; i<size; i++)  dissH2I[i] = 0.0;


  // return success
  return SUCCESS;
}

#endif
