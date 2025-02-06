/* This file implements Inflaton Potential Function */

#ifndef _MODEL_
#define _MODEL_

#include "parameters.hpp"


/* Chaotic Inflation (V = 1/2 m^2 phi^2) */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
#if MODEL==0

// Rescaling Parameters
const double rescale_B=m;
const double rescale_s=1.;


inline void modelinfo(FILE *info_)
{
  fprintf(info_,"One Field M Model\n");
  fprintf(info_,"V = 1/2 m^2 phi^2\n\n");
  fprintf(info_,"m = %e\n",m);
}

// V(phi)
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


// dV/dphi
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

// Used in Hubble Initial Condition
inline double hubble_pote()
{
  return (0.5*pw2(initfield[0]));
}


// d^2V/dphi^2
inline double dvdf2()
{
  return (1.);
}



/* Step potential Inflation (V = 1/2 m^2 phi^2 (1 + s tanh((phi-step)/d))) */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
#elif MODEL==1

// Model Parameters
const double model_s =0.01;
const double model_d =0.005;
const double model_step=14.35;

// Rescaling Parameters
const double rescale_B=m;
const double rescale_s=1.;


inline void modelinfo(FILE *info_)
{
  fprintf(info_,"One Field step potential Model\n");
  fprintf(info_,"V = 1/2 m^2 phi^2 (1 + s tanh((phi-step)/d))\n\n");
  fprintf(info_,"m = %e\n",m);
}

// V(phi)
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


// dV/dphi
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


// Used in Hubble Initial Condition
inline double hubble_pote()
{
  return (0.5*pw2(initfield[0])*(1.+model_s*tanh((initfield[0]-model_step)/model_d)))
}


// d^2/dphi^2
inline double dvdf2()
{
  return (1. + model_s*tanh((initfield[0]-model_step)/model_d)+model_s*initfield[0]*(2.-initfield[0]*tanh((initfield[0]-model_step)/model_d)/model_d)/(model_d*pw2(cosh((initfield[0]-model_step)/model_d))));
}


#endif

#endif