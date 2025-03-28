/* This file implements Inflaton Potential Function */

#ifndef _MODEL_
#define _MODEL_

/* Chaotic Inflation (V = 1/2 m^2 phi^2) */
/*-------------------------------------------------------------------------------*/
/*-------------------------------------------------------------------------------*/
#if MODEL==0

// Rescaling Parameters
const double rescale_B=m;
const double rescale_s=1.;


inline void modelinfo(FILE *info_)
{
  fprintf(info_,"Chaotic Inflation Model\n");
  fprintf(info_,"V = 1/2 m^2 phi^2\n");
  fprintf(info_,"m = %e\n",m);
}

// V(phi)
inline double potential_energy(int kloop = 1, double fsub1[nflds][N][N][N]=NULL, double fsub2[nflds][N][N][N]=NULL)
{
  DECLARE_INDICES
  double potential=0.;

  switch(kloop)
  {
    case 1:
    {
      #pragma omp parallel for reduction(+: potential) num_threads(num_thread) private(j,k)
      LOOP
      {
        potential += pw2(AXION);
      }
      potential /= pow(a,2.*rescale_s);
      break;
    }
    case 2:
    case 4:
    {
      #pragma omp parallel for reduction(+: potential) num_threads(num_thread) private(j,k)
      LOOP
      {
        potential += pw2(AXIONSUB1);
      }
      potential /= pow(asub1,2.*rescale_s);
      break;
    }
    case 3:
    {
      #pragma omp parallel for reduction(+: potential) num_threads(num_thread) private(j,k)
      LOOP
      {
        potential += pw2(AXIONSUB2);
      }
      potential /= pow(asub2,2.*rescale_s);
      break;
    }
    default:
    {
      printf("Error: kloop takes unknown value\n");
      break;
    }
  }

  potential /= (double)gridsize*2.;

  return (potential);
}


// dV/dphi
inline double dvdf(int i, int j, int k, int kloop = 1, double fsub1[nflds][N][N][N]=NULL, double fsub2[nflds][N][N][N]=NULL)
{
  double value;
  switch(kloop)
  {
    case 1:
    {
      value = AXION/pow(a,2.*rescale_s);
      break;
    }
    case 2:
    case 4:
    {
      value = AXIONSUB1/pow(asub1,2.*rescale_s);
      break;
    }
    case 3:
    {
      value = AXIONSUB2/pow(asub2,2.*rescale_s);
      break;
    }
    default:
    {
      printf("Error: kloop takes unknown value\n");
      break;
    }
  }

  return value;
}


// Used in Hubble Initial Condition
inline double hubble_pote()
{
  return (0.5*pw2(initfield[0]));
}

//d^2V/dphi^2
inline double dvdf2()
{
  return (1.);
}



/* Bumpy Axion Inflation (V = 1/2 m^2 phi^2 + Lambda^4 phi/f sin(phi/f)) */
/*-------------------------------------------------------------------------------*/
/*-------------------------------------------------------------------------------*/
#elif MODEL==1

// Model Parameters
const double couf = 1./3.3;
const double coef = 0.996;

// Rescaling Parameters
const double rescale_B=m;
const double rescale_s=1.;


inline void modelinfo(FILE *info_)
{
  fprintf(info_,"Bumpy Axion Inflation Model\n");
  fprintf(info_,"V = 1/2 m^2 phi^2 + Lambda^4 phi/f sin(phi/f))\n");
  fprintf(info_,"m = %e\n",m);
}

// V(phi)
inline double potential_energy(int kloop = 1, double fsub1[nflds][N][N][N]=NULL, double fsub2[nflds][N][N][N]=NULL)
{
  DECLARE_INDICES
  double potential=0.;

  switch(kloop)
  {
    case 1:
    {
      #pragma omp parallel for reduction(+: potential) num_threads(num_thread) private(j,k)
      LOOP
      {
        potential += pw2(AXION)+2.*coef*couf*AXION*sin(AXION/couf);
      }
      potential /= pow(a,2.*rescale_s);
      break;
    }
    case 2:
    case 4:
    {
      #pragma omp parallel for reduction(+: potential) num_threads(num_thread) private(j,k)
      LOOP
      {
        potential += pw2(AXIONSUB1)+2.*coef*couf*AXIONSUB1*sin(AXIONSUB1/couf);
      }
      potential /= pow(asub1,2.*rescale_s);
      break;
    }
    case 3:
    {
      #pragma omp parallel for reduction(+: potential) num_threads(num_thread) private(j,k)
      LOOP
      {
        potential += pw2(AXIONSUB2)+2.*coef*couf*AXIONSUB2*sin(AXIONSUB2/couf);
      }
      potential /= pow(asub2,2.*rescale_s);
      break;
    }
    default:
    {
      printf("Error: kloop takes unknown value\n");
      break;
    }
  }

  potential /= (double)gridsize*2.;

  return (potential);
}


// dV/dphi
inline double dvdf(int i, int j, int k, int kloop = 1, double fsub1[nflds][N][N][N]=NULL, double fsub2[nflds][N][N][N]=NULL)
{
  double value;
  switch(kloop)
  {
    case 1:
    {
      value = (AXION+coef*(couf*sin(AXION/couf)+AXION*cos(AXION/couf)))/pow(a,2.*rescale_s);
      break;
    }
    case 2:
    case 4:
    {
      value = (AXIONSUB1+coef*(couf*sin(AXIONSUB1/couf)+AXIONSUB1*cos(AXIONSUB1/couf)))/pow(asub1,2.*rescale_s);
      break;
    }
    case 3:
    {
      value = (AXIONSUB2+coef*(couf*sin(AXIONSUB2/couf)+AXIONSUB2*cos(AXIONSUB2/couf)))/pow(asub2,2.*rescale_s);
      break;
    }
    default:
    {
      printf("Error: kloop takes unknown value\n");
      break;
    }
  }

  return value;
}


// Used in Hubble Initial Condition
inline double hubble_pote()
{
  return (0.5*pw2(initfield[0])+coef*couf*initfield[0]*sin(initfield[0]/couf));
}


// d^2V/dphi^2
inline double dvdf2()
{
  double phase = initfield[0]/couf;

  return (1. + coef*(2.*cos(phase)-phase*sin(phase)));
}


#endif

#endif