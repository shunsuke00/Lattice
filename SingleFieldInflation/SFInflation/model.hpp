/*
General comments about the model.h file:
This file contains the following functions - all called externally
modelinfo(FILE *info_) outputs information about the model and model-specific parameters to a file.
modelinitialize() performs any model-specific initialization
potential_energy(int term, double *field_values) calculates the average potential energy density term by term. The variable num_potential_terms (just above this function) specifies how many separate potential terms are used in this model.
dvdf(int fld, int i, int j, int k) calculates the potential term in the equation of motion, dV/dfield, for the field fld at the lattice point (i,j,k)
effective_mass(double mass_sq[], double *field_values) calculates the square masses of the fields and puts them in the array mass_sq. The parameter beginning tells the function to use initial field values - if this parameter is zero then the field quantities will be calculated dynamically.
model_output(int flush, char *ext_) allows each model to include its own specialized output function(s). The parameter flush is set to 1 when infrequent calculations are being performed and 0 otherwise. The string ext_ gives the extension for output filenames.
*/
#ifndef _MODEL_
#define _MODEL_

#include "parameters.hpp"


/* Chaotic Inflation (V = 1/2 m^2 phi^2) */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
#if MODEL==0

/* - This section should be copied into parameters.h when using this model
// ---Adjustable parameters for ONEFLDM model--- //
const double m =0.51e-5;
const double m2=2.601e-11;
*/


// By default these are automatically set to A=1/f0, B=sqrt(cpl) f0^(-1+beta/2), R=6/(2+beta), S=3(2-beta)/(2+beta). They may be adjusted to different values, but the relationship S=2R-3 must be maintained for the program equations to remain correct.
const double rescale_B=m;
const double rescale_s=1.;


inline void modelinfo(FILE *info_)
{
  fprintf(info_,"One Field M Model\n");
  fprintf(info_,"V = 1/2 m^2 phi^2\n\n");
  fprintf(info_,"m = %e\n",m);
}


inline double potential_energy(int kloop=1, double fsub1[nflds][N][N][N]=NULL, double fsub2[nflds][N][N][N]=NULL)
{
  DECLARE_INDICES
  double potential=0.;

  switch(kloop)
  {
    case 1:
      #pragma omp parallel for reduction(+: potential) num_threads(num_thread) private(j,k)
      LOOP
      {
        potential += pw2(FIELD(0));
      }
      potential /= pow(a,2.*rescale_s);
      break;
    case 2:
    case 4:
      #pragma omp parallel for reduction(+: potential) num_threads(num_thread) private(j,k)
      LOOP
      {
        potential += pw2(FIELDSUB1(0));
      }
      potential /= pow(asub1,2.*rescale_s);
      break;
    case 3:
      #pragma omp parallel for reduction(+: potential) num_threads(num_thread) private(j,k)
      LOOP
      {
        potential += pw2(FIELDSUB2(0));
      }
      potential /= pow(asub2,2.*rescale_s);
      break;
  }

  potential /= (double)gridsize*2.;

  return (potential);
}


inline double dvdf(int fld, int i, int j, int k, int kloop=1, double fsub1[nflds][N][N][N]=NULL, double fsub2[nflds][N][N][N]=NULL)
{
  double value;
  switch(kloop)
  {
    case 1:
      value = FIELD(fld)/pow(a,2.*rescale_s);
      break;
    case 2:
    case 4:
      value = FIELDSUB1(fld)/pow(asub1,2.*rescale_s);
      break;
    case 3:
      value = FIELDSUB2(fld)/pow(asub2,2.*rescale_s);
      break;
  }

  return value;
}


inline double hubble_pote()
{
  return (0.5*pw2(initfield[0]));
}

inline double dvdf2()
{
  return (1.);
}




/* Step potential Inflation (V = 1/2 m^2 phi^2 (1 + s tanh((phi-step)/d))) */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
#elif MODEL==1

/* - This section should be copied into parameters.h when using this model
// ---Adjustable parameters for ONEFLD step model--- //
const double m =0.51e-5;
*/
const double model_s =0.01;
const double model_d =0.005;
const double model_step=14.35;

// By default these are automatically set to A=1/f0, B=sqrt(cpl) f0^(-1+beta/2), R=6/(2+beta), S=3(2-beta)/(2+beta). They may be adjusted to different values, but the relationship S=2R-3 must be maintained for the program equations to remain correct.
const double rescale_B=m;
const double rescale_s=1.;


inline void modelinfo(FILE *info_)
{
  fprintf(info_,"One Field step potential Model\n");
  fprintf(info_,"V = 1/2 m^2 phi^2 (1 + s tanh((phi-step)/d))\n\n");
  fprintf(info_,"m = %e\n",m);
}


inline double potential_energy(int kloop=1, double fsub1[nflds][N][N][N]=NULL, double fsub2[nflds][N][N][N]=NULL)
{
  DECLARE_INDICES
  double potential=0.;

  switch(kloop)
  {
    case 1:
      #pragma omp parallel for reduction(+: potential) num_threads(num_thread) private(j,k)
      LOOP
      {
        potential += pw2(FIELD(0))*(1. + model_s*tanh((FIELD(0)-model_step)/model_d));
      }
      potential /= pow(a,2.*rescale_s);
      break;
    case 2:
    case 4:
      #pragma omp parallel for reduction(+: potential) num_threads(num_thread) private(j,k)
      LOOP
      {
        potential += pw2(FIELDSUB1(0))*(1. + model_s*tanh((FIELDSUB1(0)-model_step)/model_d));
      }
      potential /= pow(asub1,2.*rescale_s);
      break;
    case 3:
      #pragma omp parallel for reduction(+: potential) num_threads(num_thread) private(j,k)
      LOOP
      {
        potential += pw2(FIELDSUB2(0))*(1. + model_s*tanh((FIELDSUB2(0)-model_step)/model_d));
      }
      potential /= pow(asub2,2.*rescale_s);
      break;
  }

  potential /= (double)gridsize*2.;

  return (potential);
}


inline double dvdf(int fld, int i, int j, int k, int kloop=1, double fsub1[nflds][N][N][N]=NULL, double fsub2[nflds][N][N][N]=NULL)
{
  double value;
  double phase;

  switch(kloop)
  {
    case 1:
      phase = (FIELD(fld)-model_step)/model_d;
      value = FIELD(fld)/pow(a,2.*rescale_s)*(1. + model_s*tanh(phase) + FIELD(fld)*model_s/(2.*model_d)/pw2(cosh(phase)));
      break;
    case 2:
    case 4:
      phase = (FIELDSUB1(fld)-model_step)/model_d;
      value = FIELDSUB1(fld)/pow(asub1,2.*rescale_s)*(1. + model_s*tanh(phase) + FIELDSUB1(fld)*model_s/(2.*model_d)/pw2(cosh(phase)));
      break;
    case 3:
      phase = (FIELDSUB2(fld)-model_step)/model_d;
      value = FIELDSUB2(fld)/pow(asub2,2.*rescale_s)*(1. + model_s*tanh(phase) + FIELDSUB2(fld)*model_s/(2.*model_d)/pw2(cosh(phase)));
      break;
  }

  return value;
}


inline double hubble_pote()
{
  return (0.5*pw2(initfield[0])*(1.+model_s*tanh((initfield[0]-model_step)/model_d)))
}

inline double dvdf2()
{
  return (1. + model_s*tanh((initfield[0]-model_step)/model_d)+model_s*initfield[0]*(2.-initfield[0]*tanh((initfield[0]-model_step)/model_d)/model_d)/(model_d*pw2(cosh((initfield[0]-model_step)/model_d))));
}


#endif

#endif