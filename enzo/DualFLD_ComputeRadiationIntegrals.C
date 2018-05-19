/*****************************************************************************
 * Copyright 2013 Daniel R. Reynolds                                         *
 *                                                                           *
 * This software is released under the terms of the "Enzo Public License"    *
 * in the accompanying LICENSE file.                                         *
 *****************************************************************************/
/***********************************************************************
   Two-Group (X-ray + UV) Flux-Limited Diffusion Implicit Problem Class
   Radiation Spectrum Evaluation routine 

   author: Daniel R. Reynolds
   date:   January 2013

   PURPOSE: Computes the UV and Xray radiation spectrum integrals
                int_{nu0_HI}^{inf} chi_E(nu) dnu
                int_{nu0_HI}^{inf} chi_E(nu)*sigma_HI(nu) dnu
                int_{nu0_HI}^{inf} chi_E(nu)*sigma_HI(nu)/nu dnu
                int_{nu0_HeI}^{inf} chi_E(nu)*sigma_HeI(nu) dnu
                int_{nu0_HeI}^{inf} chi_E(nu)*sigma_HeI(nu)/nu dnu
                int_{nu0_HeII}^{inf} chi_E(nu)*sigma_HeII(nu) dnu
                int_{nu0_HeII}^{inf} chi_E(nu)*sigma_HeII(nu)/nu dnu
            where nu0_* is the ionization threshold of the relevant 
            species, chi_E(nu) is the spectrum of the relevant radiation 
            energy density, sigma_HI(nu) is the ionization cross section 
            of HI, sigma_HeI(nu) is the ionization cross section of HeI, 
            and sigma_HeII(nu) is the ionization cross section of HeII.

            These are computed using a simple 4th-order accurate 
            numerical quadrature rule (composite Simpson's), on the 
            re-mapped indefinite integral
                   int_{nu0}^{inf} f(nu) dnu
                 = int_0^1 nu0/(x^2)*f(nu0/x) dx
            The composite Simpson's rule computes this integral as 
            the sum of small quadrature intervals [xl,xr]_i, where 
                   [0,1] = Union_i [xl,xr]_i,
            and xm = (xl+xr)/2, via the quadrature formula
                   int_{xl}^{xr} g(x) dx 
                 = (xr-xl)/6*(g(xl) + 4*g(xm) + g(xr))
            Note: since 1/0 = infty, we begin the mapped quadrature 
            integration at a small positive value.

            We further note that these integrals are typically 
            computed once and stored in the DualFLD object, upon 
            initialization of said object.  We compute these here to 
            allow for control over the assumed radiation energy 
            density spectrum, as stored in 
            DualFLD_RadiationSpectrum.C, as opposed to fixing the 
            spectrum, calculating these same integrals outside of the 
            simulation, and hard-coding the physics routines to use 
            the fixed integral values.
************************************************************************/
#ifdef RAD_HYDRO
#include "DualFLD.h"


int DualFLD::ComputeRadiationIntegrals() {

  if (debug)  printf("Entering DualFLD::ComputeRadiationIntegrals\n");

  // set necessary constants
  float h        = 6.6260693e-27;        // Planck's constant [ergs*s]
  float ev2erg   = 1.60217653e-12;       // conversion constant from eV to ergs
  float nu0_HI   = hnu0_HI*ev2erg/h;     // ionization threshold of HI (hz)
  float nu0_HeI  = hnu0_HeI*ev2erg/h;    // ionization threshold of HeI (hz)
  float nu0_HeII = hnu0_HeII*ev2erg/h;   // ionization threshold of HeII (hz)
  float epsilon  = 1.0;                  // floating point roundoff
  while ((1.0 + epsilon) > 1.0)  epsilon*=0.5;

  // set integration parameters
  float FreqH = 1e-4;             // normalized width of quadrature bins
  float Llimit = POW(epsilon,0.35);     // integration limits, shift away from 
  float Ulimit = 1.0-POW(epsilon,0.47); //   [0,1] to avoid over/underflow errors

  // initialize quadrature points and values
  float xl, xm, xr, nu_l, nu_m, nu_r, fl, fm, fr;
  float fl_E, fm_E, fr_E, fl_ni, fm_ni, fr_ni, fl_nu, fm_nu, fr_nu;
  int i;


  ////////////////////////////////////////////////////////////////
  // Xray spectrum calculations
  {
    for (i=0; i<7; i++)  RadIntXr[i] = 0.0;

    // integrate based on XrSpectrum parameter, since for 
    // monochromatic problems integration uses delta function 
    switch (XrSpectrum) {

    case -1:     // monochromatic at XrFrequency

      // evaluation point
      nu_m = XrFrequency*ev2erg/h*(1.0 + 2.0*epsilon);

      // integral( delta_{nu0}(nu) )
      RadIntXr[0] = 1.0;

      // integral( delta_{nu0}(nu) * sig_{HI}(nu) )
      if (XrFrequency >= hnu0_HI)
	RadIntXr[1] = this->CrossSections(nu_m,0);

      // integral( delta_{nu0}(nu) * sig_{HI}(nu) / nu )
      if (XrFrequency >= hnu0_HI)
	RadIntXr[2] = this->CrossSections(nu_m,0) / nu_m;

      // integral( delta_{nu0}(nu) * sig_{HeI}(nu) ) 
      if (XrFrequency >= hnu0_HeI)
	RadIntXr[3] = this->CrossSections(nu_m,1);

      // integral( delta_{nu0}(nu) * sig_{HeI}(nu) / nu )
      if (XrFrequency >= hnu0_HeI)
	RadIntXr[4] = this->CrossSections(nu_m,1) / nu_m;

      // integral( delta_{nu0}(nu) * sig_{HeII}(nu) )
      if (XrFrequency >= hnu0_HeII)
	RadIntXr[5] = this->CrossSections(nu_m,2);

      // integral( delta_{nu0}(nu) * sig_{HeII}(nu) / nu )
      if (XrFrequency >= hnu0_HeII)
	RadIntXr[6] = this->CrossSections(nu_m,2) / nu_m;

      break;


    default:     // numerically integrate 

      //////////
      // Compute intChiE, intChiESigHI and intChiESigHInu integrals

      //   left end of integration
      //      can't start at 0, shift over a bit
      xr   = Llimit;
      nu_r = nu0_HI/xr;
      //      get function values at this location
      fr_E  = this->RadiationSpectrum(XrSpectrum, nu_r);
      fr_ni = this->CrossSections(nu_r,0);
      fr_nu = 1.0/nu_r;
      if ((fr_E == -1.0) || (fr_ni == -1.0)) {
	fprintf(stderr,"ComputeRadiationIntegrals Error in evaluating spectrum\n");
	return FAIL;
      }
  
      //   iterate over intervals
      for (i=1; i<1e9; i++) {

	//      set quadrature points in interval
	xl = xr;  // cannot start at 0, so shift over a bit
	xr = min(xl+FreqH,Ulimit);
	xm = 0.5*(xl+xr);

	//      copy left subinterval function value, location, etc
	nu_l  = nu_r;
	fl_E  = fr_E;
	fl_ni = fr_ni;
	fl_nu = fr_nu;

	//      evaluate chi_E(), sigma_HI, 1/nu at remapped quad. pts.
	nu_m  = nu0_HI/xm;
	fm_E  = this->RadiationSpectrum(XrSpectrum, nu_m);
	fm_ni = this->CrossSections(nu_m,0);
	fm_nu = 1.0/nu_m;
	if ((fm_E == -1.0) || (fm_ni == -1.0)) {
	  fprintf(stderr,"ComputeRadiationIntegrals Error in evaluating spectrum\n");
	  return FAIL;
	}
  
	nu_r  = nu0_HI/xr;
	fr_E  = this->RadiationSpectrum(XrSpectrum, nu_r);
	fr_ni = this->CrossSections(nu_r,0);
	fr_nu = 1.0/nu_r;
	if ((fr_E == -1.0) || (fr_ni == -1.0)) {
	  fprintf(stderr,"ComputeRadiationIntegrals Error in evaluating spectrum\n");
	  return FAIL;
	}
  
	//      compute integrals using re-scaled function values
	fl = fl_E/xl/xl;
	fm = fm_E/xm/xm;
	fr = fr_E/xr/xr;
	RadIntXr[0] += nu0_HI*(xr-xl)/6.0*(fl + 4.0*fm + fr);
    
	fl = fl_E*fl_ni/xl/xl;
	fm = fm_E*fm_ni/xm/xm;
	fr = fr_E*fr_ni/xr/xr;  
	RadIntXr[1] += nu0_HI*(xr-xl)/6.0*(fl + 4.0*fm + fr);
    
	fl = fl_E*fl_ni*fl_nu/xl/xl;
	fm = fm_E*fm_ni*fm_nu/xm/xm;
	fr = fr_E*fr_ni*fr_nu/xr/xr;
	RadIntXr[2] += nu0_HI*(xr-xl)/6.0*(fl + 4.0*fm + fr);

	//      quit if we have finished interval
	if (xr >= (Ulimit-10.0*epsilon))  break;
  
      }
  

      //////////
      // Compute intChiESigHeI and intChiESigHeInu integrals

      //   left end of integration
      //      can't start at 0, shift over a bit
      xr   = Llimit;
      nu_r = nu0_HeI/xr;
      //      get function values at this location
      fr_E  = this->RadiationSpectrum(XrSpectrum, nu_r);
      fr_ni = this->CrossSections(nu_r,1);
      fr_nu = 1.0/nu_r;
      if ((fr_E == -1.0) || (fr_ni == -1.0)) {
	fprintf(stderr,"ComputeRadiationIntegrals Error in evaluating spectrum\n");
	return FAIL;
      }
  
      //   iterate over intervals
      for (i=1; i<1e9; i++) {

	//      set quadrature points in interval
	xl = xr;  // cannot start at 0, so shift over a bit
	xr = min(xl+FreqH,Ulimit);
	xm = 0.5*(xl+xr);

	//      copy left subinterval function value, location, etc
	nu_l  = nu_r;
	fl_E  = fr_E;
	fl_ni = fr_ni;
	fl_nu = fr_nu;

	//      evaluate chi_E(), sigma_HeI, 1/nu at remapped quad. pts.
	nu_m  = nu0_HeI/xm;
	fm_E  = this->RadiationSpectrum(XrSpectrum, nu_m);
	fm_ni = this->CrossSections(nu_m,1);
	fm_nu = 1.0/nu_m;
	if ((fm_E == -1.0) || (fm_ni == -1.0)) {
	  fprintf(stderr,"ComputeRadiationIntegrals Error in evaluating spectrum\n");
	  return FAIL;
	}
  
	nu_r  = nu0_HeI/xr;
	fr_E  = this->RadiationSpectrum(XrSpectrum, nu_r);
	fr_ni = this->CrossSections(nu_r,1);
	fr_nu = 1.0/nu_r;
	if ((fr_E == -1.0) || (fr_ni == -1.0)) {
	  fprintf(stderr,"ComputeRadiationIntegrals Error in evaluating spectrum\n");
	  return FAIL;
	}
  
	//      compute integrals using re-scaled function values
	fl = fl_E*fl_ni/xl/xl;
	fm = fm_E*fm_ni/xm/xm;
	fr = fr_E*fr_ni/xr/xr;
	RadIntXr[3] += nu0_HeI*(xr-xl)/6.0*(fl + 4.0*fm + fr);
    
	fl = fl_E*fl_ni*fl_nu/xl/xl;
	fm = fm_E*fm_ni*fm_nu/xm/xm;
	fr = fr_E*fr_ni*fr_nu/xr/xr;
	RadIntXr[4] += nu0_HeI*(xr-xl)/6.0*(fl + 4.0*fm + fr);

	//      quit if we have finished interval
	if (xr >= (Ulimit-10.0*epsilon))  break;
  
      }
  

      //////////
      // Compute intChiESigHeII and intChiESigHeIInu integrals

      //   left end of integration
      //      can't start at 0, shift over a bit
      xr   = Llimit;
      nu_r = nu0_HeII/xr;
      //      get function values at this location
      fr_E  = this->RadiationSpectrum(XrSpectrum, nu_r);
      fr_ni = this->CrossSections(nu_r,2);
      fr_nu = 1.0/nu_r;
      if ((fr_E == -1.0) || (fr_ni == -1.0)) {
	fprintf(stderr,"ComputeRadiationIntegrals Error in evaluating spectrum\n");
	return FAIL;
      }
  
      //   iterate over intervals
      for (i=1; i<1e9; i++) {

	//      set quadrature points in interval
	xl = xr;  // cannot start at 0, so shift over a bit
	xr = min(xl+FreqH,Ulimit);
	xm = 0.5*(xl+xr);

	//      copy left subinterval function value, location, etc
	nu_l  = nu_r;
	fl_E  = fr_E;
	fl_ni = fr_ni;
	fl_nu = fr_nu;

	//      evaluate chi_E(), sigma_HeII, 1/nu at remapped quad. pts.
	nu_m  = nu0_HeII/xm;
	fm_E  = this->RadiationSpectrum(XrSpectrum, nu_m);
	fm_ni = this->CrossSections(nu_m,2);
	fm_nu = 1.0/nu_m;
	if ((fm_E == -1.0) || (fm_ni == -1.0)) {
	  fprintf(stderr,"ComputeRadiationIntegrals Error in evaluating spectrum\n");
	  return FAIL;
	}
  
	nu_r  = nu0_HeII/xr;
	fr_E  = this->RadiationSpectrum(XrSpectrum, nu_r);
	fr_ni = this->CrossSections(nu_r,2);
	fr_nu = 1.0/nu_r;
	if ((fr_E == -1.0) || (fr_ni == -1.0)) {
	  fprintf(stderr,"ComputeRadiationIntegrals Error in evaluating spectrum\n");
	  return FAIL;
	}
  
	//      compute integrals using re-scaled function values
	fl = fl_E*fl_ni/xl/xl;
	fm = fm_E*fm_ni/xm/xm;
	fr = fr_E*fr_ni/xr/xr;
	RadIntXr[5] += nu0_HeII*(xr-xl)/6.0*(fl + 4.0*fm + fr);
    
	fl = fl_E*fl_ni*fl_nu/xl/xl;
	fm = fm_E*fm_ni*fm_nu/xm/xm;
	fr = fr_E*fr_ni*fr_nu/xr/xr;
	RadIntXr[6] += nu0_HeII*(xr-xl)/6.0*(fl + 4.0*fm + fr);

	//      quit if we have finished interval
	if (xr >= (Ulimit-10.0*epsilon))  break;
  
      }

    }  // end switch(XrSpectrum)

    if (debug) {
      printf("  Computed X-ray Radiation Integrals:\n");
      printf("    intChiE          = %22.16e\n",RadIntXr[0]);
      printf("    intChiESigHI     = %22.16e\n",RadIntXr[1]);
      printf("    intChiESigHInu   = %22.16e\n",RadIntXr[2]);
      printf("    intChiESigHeI    = %22.16e\n",RadIntXr[3]);
      printf("    intChiESigHeInu  = %22.16e\n",RadIntXr[4]);
      printf("    intChiESigHeII   = %22.16e\n",RadIntXr[5]);
      printf("    intChiESigHeIInu = %22.16e\n",RadIntXr[6]);
    }

  }  // end Xray spectrum calculations





  ////////////////////////////////////////////////////////////////
  // UV spectrum calculations
  {
    for (i=0; i<7; i++)  RadIntUV[i] = 0.0;

    // integrate based on UVSpectrum parameter, since for 
    // monochromatic problems integration uses delta function 
    switch (UVSpectrum) {

    case -1:      // monochromatic at UVFrequency

      // evaluation point
      nu_m = UVFrequency*ev2erg/h*(1.0 + 2.0*epsilon);
      
      // integral( delta_{nu0}(nu) )
      RadIntUV[0] = 1.0;
      
      // integral( delta_{nu0}(nu) * sig_{HI}(nu) )
      if (UVFrequency >= hnu0_HI)
	RadIntUV[1] = this->CrossSections(nu_m,0);
      
      // integral( delta_{nu0}(nu) * sig_{HI}(nu) / nu )
      if (UVFrequency >= hnu0_HI)
	RadIntUV[2] = this->CrossSections(nu_m,0) / nu_m;

      // integral( delta_{nu0}(nu) * sig_{HeI}(nu) ) 
      if (UVFrequency >= hnu0_HeI)
	RadIntUV[3] = this->CrossSections(nu_m,1);

      // integral( delta_{nu0}(nu) * sig_{HeI}(nu) / nu )
      if (UVFrequency >= hnu0_HeI)
	RadIntUV[4] = this->CrossSections(nu_m,1) / nu_m;

      // integral( delta_{nu0}(nu) * sig_{HeII}(nu) )
      if (UVFrequency >= hnu0_HeII)
	RadIntUV[5] = this->CrossSections(nu_m,2);

      // integral( delta_{nu0}(nu) * sig_{HeII}(nu) / nu )
      if (UVFrequency >= hnu0_HeII)
	RadIntUV[6] = this->CrossSections(nu_m,2) / nu_m;

      break;

      
    default:     // numerically integrate 

      //////////
      // Compute intChiE, intChiESigHI and intChiESigHInu integrals

      //   left end of integration
      //      can't start at 0, shift over a bit
      xr   = Llimit;
      nu_r = nu0_HI/xr;
      //      get function values at this location
      fr_E  = this->RadiationSpectrum(UVSpectrum, nu_r);
      fr_ni = this->CrossSections(nu_r,0);
      fr_nu = 1.0/nu_r;
      if ((fr_E == -1.0) || (fr_ni == -1.0)) {
	fprintf(stderr,"ComputeRadiationIntegrals Error in evaluating spectrum\n");
	return FAIL;
      }
  
      //   iterate over intervals
      for (i=1; i<1e9; i++) {

	//      set quadrature points in interval
	xl = xr;  // cannot start at 0, so shift over a bit
	xr = min(xl+FreqH,Ulimit);
	xm = 0.5*(xl+xr);

	//      copy left subinterval function value, location, etc
	nu_l  = nu_r;
	fl_E  = fr_E;
	fl_ni = fr_ni;
	fl_nu = fr_nu;

	//      evaluate chi_E(), sigma_HI, 1/nu at remapped quad. pts.
	nu_m  = nu0_HI/xm;
	fm_E  = this->RadiationSpectrum(UVSpectrum, nu_m);
	fm_ni = this->CrossSections(nu_m,0);
	fm_nu = 1.0/nu_m;
	if ((fm_E == -1.0) || (fm_ni == -1.0)) {
	  fprintf(stderr,"ComputeRadiationIntegrals Error in evaluating spectrum\n");
	  return FAIL;
	}
  
	nu_r  = nu0_HI/xr;
	fr_E  = this->RadiationSpectrum(UVSpectrum, nu_r);
	fr_ni = this->CrossSections(nu_r,0);
	fr_nu = 1.0/nu_r;
	if ((fr_E == -1.0) || (fr_ni == -1.0)) {
	  fprintf(stderr,"ComputeRadiationIntegrals Error in evaluating spectrum\n");
	  return FAIL;
	}
  
	//      compute integrals using re-scaled function values
	fl = fl_E/xl/xl;
	fm = fm_E/xm/xm;
	fr = fr_E/xr/xr;
	RadIntUV[0] += nu0_HI*(xr-xl)/6.0*(fl + 4.0*fm + fr);
    
	fl = fl_E*fl_ni/xl/xl;
	fm = fm_E*fm_ni/xm/xm;
	fr = fr_E*fr_ni/xr/xr;  
	RadIntUV[1] += nu0_HI*(xr-xl)/6.0*(fl + 4.0*fm + fr);
    
	fl = fl_E*fl_ni*fl_nu/xl/xl;
	fm = fm_E*fm_ni*fm_nu/xm/xm;
	fr = fr_E*fr_ni*fr_nu/xr/xr;
	RadIntUV[2] += nu0_HI*(xr-xl)/6.0*(fl + 4.0*fm + fr);

	//      quit if we have finished interval
	if (xr >= (Ulimit-10.0*epsilon))  break;
  
      }
  

      //////////
      // Compute intChiESigHeI and intChiESigHeInu integrals

      //   left end of integration
      //      can't start at 0, shift over a bit
      xr   = Llimit;
      nu_r = nu0_HeI/xr;
      //      get function values at this location
      fr_E  = this->RadiationSpectrum(UVSpectrum, nu_r);
      fr_ni = this->CrossSections(nu_r,1);
      fr_nu = 1.0/nu_r;
      if ((fr_E == -1.0) || (fr_ni == -1.0)) {
	fprintf(stderr,"ComputeRadiationIntegrals Error in evaluating spectrum\n");
	return FAIL;
      }
  
      //   iterate over intervals
      for (i=1; i<1e9; i++) {

	//      set quadrature points in interval
	xl = xr;  // cannot start at 0, so shift over a bit
	xr = min(xl+FreqH,Ulimit);
	xm = 0.5*(xl+xr);

	//      copy left subinterval function value, location, etc
	nu_l  = nu_r;
	fl_E  = fr_E;
	fl_ni = fr_ni;
	fl_nu = fr_nu;

	//      evaluate chi_E(), sigma_HeI, 1/nu at remapped quad. pts.
	nu_m  = nu0_HeI/xm;
	fm_E  = this->RadiationSpectrum(UVSpectrum, nu_m);
	fm_ni = this->CrossSections(nu_m,1);
	fm_nu = 1.0/nu_m;
	if ((fm_E == -1.0) || (fm_ni == -1.0)) {
	  fprintf(stderr,"ComputeRadiationIntegrals Error in evaluating spectrum\n");
	  return FAIL;
	}
  
	nu_r  = nu0_HeI/xr;
	fr_E  = this->RadiationSpectrum(UVSpectrum, nu_r);
	fr_ni = this->CrossSections(nu_r,1);
	fr_nu = 1.0/nu_r;
	if ((fr_E == -1.0) || (fr_ni == -1.0)) {
	  fprintf(stderr,"ComputeRadiationIntegrals Error in evaluating spectrum\n");
	  return FAIL;
	}
  
	//      compute integrals using re-scaled function values
	fl = fl_E*fl_ni/xl/xl;
	fm = fm_E*fm_ni/xm/xm;
	fr = fr_E*fr_ni/xr/xr;
	RadIntUV[3] += nu0_HeI*(xr-xl)/6.0*(fl + 4.0*fm + fr);
    
	fl = fl_E*fl_ni*fl_nu/xl/xl;
	fm = fm_E*fm_ni*fm_nu/xm/xm;
	fr = fr_E*fr_ni*fr_nu/xr/xr;
	RadIntUV[4] += nu0_HeI*(xr-xl)/6.0*(fl + 4.0*fm + fr);

	//      quit if we have finished interval
	if (xr >= (Ulimit-10.0*epsilon))  break;
  
      }
  

      //////////
      // Compute intChiESigHeII and intChiESigHeIInu integrals

      //   left end of integration
      //      can't start at 0, shift over a bit
      xr   = Llimit;
      nu_r = nu0_HeII/xr;
      //      get function values at this location
      fr_E  = this->RadiationSpectrum(UVSpectrum, nu_r);
      fr_ni = this->CrossSections(nu_r,2);
      fr_nu = 1.0/nu_r;
      if ((fr_E == -1.0) || (fr_ni == -1.0)) {
	fprintf(stderr,"ComputeRadiationIntegrals Error in evaluating spectrum\n");
	return FAIL;
      }
  
      //   iterate over intervals
      for (i=1; i<1e9; i++) {

	//      set quadrature points in interval
	xl = xr;  // cannot start at 0, so shift over a bit
	xr = min(xl+FreqH,Ulimit);
	xm = 0.5*(xl+xr);

	//      copy left subinterval function value, location, etc
	nu_l  = nu_r;
	fl_E  = fr_E;
	fl_ni = fr_ni;
	fl_nu = fr_nu;

	//      evaluate chi_E(), sigma_HeII, 1/nu at remapped quad. pts.
	nu_m  = nu0_HeII/xm;
	fm_E  = this->RadiationSpectrum(UVSpectrum, nu_m);
	fm_ni = this->CrossSections(nu_m,2);
	fm_nu = 1.0/nu_m;
	if ((fm_E == -1.0) || (fm_ni == -1.0)) {
	  fprintf(stderr,"ComputeRadiationIntegrals Error in evaluating spectrum\n");
	  return FAIL;
	}
  
	nu_r  = nu0_HeII/xr;
	fr_E  = this->RadiationSpectrum(UVSpectrum, nu_r);
	fr_ni = this->CrossSections(nu_r,2);
	fr_nu = 1.0/nu_r;
	if ((fr_E == -1.0) || (fr_ni == -1.0)) {
	  fprintf(stderr,"ComputeRadiationIntegrals Error in evaluating spectrum\n");
	  return FAIL;
	}
  
	//      compute integrals using re-scaled function values
	fl = fl_E*fl_ni/xl/xl;
	fm = fm_E*fm_ni/xm/xm;
	fr = fr_E*fr_ni/xr/xr;
	RadIntUV[5] += nu0_HeII*(xr-xl)/6.0*(fl + 4.0*fm + fr);
    
	fl = fl_E*fl_ni*fl_nu/xl/xl;
	fm = fm_E*fm_ni*fm_nu/xm/xm;
	fr = fr_E*fr_ni*fr_nu/xr/xr;
	RadIntUV[6] += nu0_HeII*(xr-xl)/6.0*(fl + 4.0*fm + fr);

	//      quit if we have finished interval
	if (xr >= (Ulimit-10.0*epsilon))  break;
  
      }

    }  // end switch(UVSpectrum)

    if (debug) {
      printf("  Computed UV Radiation Integrals:\n");
      printf("    intChiE          = %22.16e\n",RadIntUV[0]);
      printf("    intChiESigHI     = %22.16e\n",RadIntUV[1]);
      printf("    intChiESigHInu   = %22.16e\n",RadIntUV[2]);
      printf("    intChiESigHeI    = %22.16e\n",RadIntUV[3]);
      printf("    intChiESigHeInu  = %22.16e\n",RadIntUV[4]);
      printf("    intChiESigHeII   = %22.16e\n",RadIntUV[5]);
      printf("    intChiESigHeIInu = %22.16e\n",RadIntUV[6]);
    }

  }  // end UV spectrum calculations


  return SUCCESS;
}

#endif
