#ifndef _ELEMENT_EVO_
#define _ELEMENT_EVO_

#include "latticeeasy.hpp"

/* Heavily used functions */
/*------------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------------*/
inline int INCREMENT(int i)
{
  return( (i==N-1) ? 0 : i+1 );
}

inline int DECREMENT(int i)
{
  return( (i==0) ? N-1 : i-1 );
}


/* Laplacian */
/*------------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------------*/
inline double lapl(int fld, int i, int j, int k, int kloop_arg, ARG)
{
  double value;

  double (*f)[N][N][N] = flist[0];
  double (*fsub1)[N][N][N] = flist[2];
  double (*fsub2)[N][N][N] = flist[4];

  switch(kloop_arg)
  {
    case 1:
      value = (f[fld][i][j][k+1] + f[fld][i][j][k-1]
            +f[fld][i][j+1][k] + f[fld][i][j-1][k]
            +f[fld][i+1][j][k] + f[fld][i-1][j][k]
            -6.*f[fld][i][j][k]);
      break;
    case 2:
    case 4:
      value = (fsub1[fld][i][j][k+1] + fsub1[fld][i][j][k-1]
            +fsub1[fld][i][j+1][k] + fsub1[fld][i][j-1][k]
            +fsub1[fld][i+1][j][k] + fsub1[fld][i-1][j][k]
            -6.*fsub1[fld][i][j][k]);
      break;
    case 3:
      value = (fsub2[fld][i][j][k+1] + fsub2[fld][i][j][k-1]
            +fsub2[fld][i][j+1][k] + fsub2[fld][i][j-1][k]
            +fsub2[fld][i+1][j][k] + fsub2[fld][i-1][j][k]
            -6.*fsub2[fld][i][j][k]);
      break;
  }

  return value;
}


inline double laplb(int fld, int i, int j, int k, int kloop_arg, ARG)
{
  double value;

  double (*f)[N][N][N] = flist[0];
  double (*fsub1)[N][N][N] = flist[2];
  double (*fsub2)[N][N][N] = flist[4];

  switch(kloop_arg)
  {
    case 1:
      value = (f[fld][i][j][INCREMENT(k)] + f[fld][i][j][DECREMENT(k)]
            +f[fld][i][INCREMENT(j)][k] + f[fld][i][DECREMENT(j)][k]
            +f[fld][INCREMENT(i)][j][k] + f[fld][DECREMENT(i)][j][k]
            -6.*f[fld][i][j][k]);
      break;
    case 2:
    case 4:
      value = (fsub1[fld][i][j][INCREMENT(k)] + fsub1[fld][i][j][DECREMENT(k)]
            +fsub1[fld][i][INCREMENT(j)][k] + fsub1[fld][i][DECREMENT(j)][k]
            +fsub1[fld][INCREMENT(i)][j][k] + fsub1[fld][DECREMENT(i)][j][k]
            -6.*fsub1[fld][i][j][k]);
      break;
    case 3:
      value = (fsub2[fld][i][j][INCREMENT(k)] + fsub2[fld][i][j][DECREMENT(k)]
            +fsub2[fld][i][INCREMENT(j)][k] + fsub2[fld][i][DECREMENT(j)][k]
            +fsub2[fld][INCREMENT(i)][j][k] + fsub2[fld][DECREMENT(i)][j][k]
            -6.*fsub2[fld][i][j][k]);
      break;
  }

  return value;
}


/* Gradient energy (externally called) 1/2 <|Grad(f)|^2> = 1/2 <-f Lapl(f)> */
/*------------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------------*/
double gradient_energy(int fld, int kloop_arg, ARG)
{
  int i=0,j=0,k=0;
  double gradient=0.;
  double norm=1./pw2(dx); // Converts lapl() to actual Laplacian

  double (*f)[N][N][N] = flist[0];
  double (*fsub1)[N][N][N] = flist[2];
  double (*fsub2)[N][N][N] = flist[4];

  switch(kloop_arg)
  {
    case 1:
    #pragma omp parallel for reduction(-:gradient) num_threads(num_thread) private(j,k)
      for(i=1;i<N-1;i++)
        for(j=1;j<N-1;j++)
          for(k=1;k<N-1;k++)
          {
            gradient -= FIELD(fld)*lapl(fld,i,j,k,kloop_arg,flist);
          }
    #pragma omp parallel for reduction(-:gradient) num_threads(num_thread) private(j)
      for(i=0;i<N;i++)
        for(j=0;j<N;j++)
        {
          gradient -= f[fld][i][j][0]*laplb(fld,i,j,0,kloop_arg,flist) + f[fld][i][j][N-1]*laplb(fld,i,j,N-1,kloop_arg,flist); // z=0,N-1
          if(j==0 || j==N-1) continue;
          gradient -= f[fld][i][0][j]*laplb(fld,i,0,j,kloop_arg,flist) + f[fld][i][N-1][j]*laplb(fld,i,N-1,j,kloop_arg,flist); // y=0,N-1
          if(i==0 || i==N-1) continue;
          gradient -= f[fld][0][i][j]*laplb(fld,0,i,j,kloop_arg,flist) + f[fld][N-1][i][j]*laplb(fld,N-1,i,j,kloop_arg,flist); // x=0,N-1
        }
      break;
    case 2:
    case 4:
    #pragma omp parallel for reduction(-:gradient) num_threads(num_thread) private(j,k)
      for(i=1;i<N-1;i++)
        for(j=1;j<N-1;j++)
          for(k=1;k<N-1;k++)
          {
            gradient -= FIELDSUB1(fld)*lapl(fld,i,j,k,kloop_arg,flist);
          }
    #pragma omp parallel for reduction(-:gradient) num_threads(num_thread) private(j)
      for(i=0;i<N;i++)
        for(j=0;j<N;j++)
        {
          gradient -= fsub1[fld][i][j][0]*laplb(fld,i,j,0,kloop_arg,flist) + fsub1[fld][i][j][N-1]*laplb(fld,i,j,N-1,kloop_arg,flist); // z=0,N-1
          if(j==0 || j==N-1) continue;
          gradient -= fsub1[fld][i][0][j]*laplb(fld,i,0,j,kloop_arg,flist) + fsub1[fld][i][N-1][j]*laplb(fld,i,N-1,j,kloop_arg,flist); // y=0,N-1
          if(i==0 || i==N-1) continue;
          gradient -= fsub1[fld][0][i][j]*laplb(fld,0,i,j,kloop_arg,flist) + fsub1[fld][N-1][i][j]*laplb(fld,N-1,i,j,kloop_arg,flist); // x=0,N-1
        }
      break;
    case 3:
    #pragma omp parallel for reduction(-:gradient) num_threads(num_thread) private(j,k)
      for(i=1;i<N-1;i++)
        for(j=1;j<N-1;j++)
          for(k=1;k<N-1;k++)
          {
            gradient -= FIELDSUB2(fld)*lapl(fld,i,j,k,kloop_arg,flist);
          }
    #pragma omp parallel for reduction(-:gradient) num_threads(num_thread) private(j)    
      for(i=0;i<N;i++)
        for(j=0;j<N;j++)
        {
          gradient -= fsub2[fld][i][j][0]*laplb(fld,i,j,0,kloop_arg,flist) + fsub2[fld][i][j][N-1]*laplb(fld,i,j,N-1,kloop_arg,flist); // z=0,N-1
          if(j==0 || j==N-1) continue;
          gradient -= fsub2[fld][i][0][j]*laplb(fld,i,0,j,kloop_arg,flist) + fsub2[fld][i][N-1][j]*laplb(fld,i,N-1,j,kloop_arg,flist); // y=0,N-1
          if(i==0 || i==N-1) continue;
          gradient -= fsub2[fld][0][i][j]*laplb(fld,0,i,j,kloop_arg,flist) + fsub2[fld][N-1][i][j]*laplb(fld,N-1,i,j,kloop_arg,flist); // x=0,N-1
        }
  }

  return(.5*gradient*norm/(double)gridsize);
}


/* rho-3p (used in ad evolution) */
/*------------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------------*/
double stress_3pressure(int kloop_arg, ARG)
{
  int i,j,k,fld;
  double rho_3p=0.;

  double (*fd)[N][N][N] = flist[1];
  double (*fdsub1)[N][N][N] = flist[3];
  double (*fdsub2)[N][N][N] = flist[5];

  switch(kloop_arg)
  {
    case 1:
    #pragma omp parallel for reduction(-:rho_3p) num_threads(num_thread) private(j,k)
      LOOP
      {
        rho_3p -= pw2(FIELDD(0));
      }
      rho_3p /= (double)gridsize*pw2(a);
      rho_3p += 2.*gradient_energy(0,kloop_arg,flist)/pow(a,2.*rescale_s+2.);
      rho_3p += 4.*potential_energy(kloop_arg,flist);
      break;
    case 2:
    case 4:
    #pragma omp parallel for reduction(-:rho_3p) num_threads(num_thread) private(j,k)
      LOOP
      {
        rho_3p -= pw2(FIELDDSUB1(0));
      }
      rho_3p /= (double)gridsize*pw2(asub1);
      rho_3p += 2.*gradient_energy(0,kloop_arg,flist)/pow(asub1,2.*rescale_s+2.);
      rho_3p += 4.*potential_energy(kloop_arg,flist);
      break;
    case 3:
    #pragma omp parallel for reduction(-:rho_3p) num_threads(num_thread) private(j,k)
      LOOP
      {
        rho_3p -= pw2(FIELDDSUB2(0));
      }
      rho_3p /= (double)gridsize*pw2(asub2);
      rho_3p += 2.*gradient_energy(0,kloop_arg,flist)/pow(asub2,2.*rescale_s+2.);
      rho_3p += 4.*potential_energy(kloop_arg,flist);
      break;
  }

  return rho_3p;
}



#endif