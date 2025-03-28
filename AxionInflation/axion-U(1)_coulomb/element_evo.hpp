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

inline int INCREMENT2(int i)
{
  return( (i>=N-2) ? ( (i==N-1) ? 1 : 0 ) : i+2 );
}

inline int DECREMENT2(int i)
{
  return( (i<=1) ? ( (i==0) ? N-2 : N-1 ) : i-2 );
}


/* Single Derivative */
/*------------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------------*/
// x-Direction
inline double deriv_x(int fld, int i, int j, int k, int kloop_arg, FSUB_ARG)
{
  double value;

  switch(kloop_arg)
  {
    case 1:
    {
      value = (f[fld][i+1][j][k] - f[fld][i-1][j][k]);
      break;
    }
    case 2:
    case 4:
    {
      value = (fsub1[fld][i+1][j][k] - fsub1[fld][i-1][j][k]);
      break;
    }
    case 3:
    {
      value = (fsub2[fld][i+1][j][k] - fsub2[fld][i-1][j][k]);
      break;
    }
    default:
    {
      printf("Error: kloop takes unknown value\n");
      break;
    }
  }
  value = value/(2.*dx);

  return value;
}

inline double deriv_xb(int fld, int i, int j, int k, int kloop_arg, FSUB_ARG)
{
  double value;

  switch(kloop_arg)
  {
    case 1:
    {
      value = (f[fld][INCREMENT(i)][j][k] - f[fld][DECREMENT(i)][j][k]);
      break;
    }
    case 2:
    case 4:
    {
      value = (fsub1[fld][INCREMENT(i)][j][k] - fsub1[fld][DECREMENT(i)][j][k]);
      break;
    }
    case 3:
    {
      value = (fsub2[fld][INCREMENT(i)][j][k] - fsub2[fld][DECREMENT(i)][j][k]);
      break;
    }
    default:
    {
      printf("Error: kloop takes unknown value\n");
      break;
    }
  }
  value = value/(2.*dx);

  return value;
}

// y-Direction
inline double deriv_y(int fld, int i, int j, int k, int kloop_arg, FSUB_ARG)
{
  double value;

  switch(kloop_arg)
  {
    case 1:
    {
      value = (f[fld][i][j+1][k] - f[fld][i][j-1][k]);
      break;
    }
    case 2:
    case 4:
    {
      value = (fsub1[fld][i][j+1][k] - fsub1[fld][i][j-1][k]);
      break;
    }
    case 3:
    {
      value = (fsub2[fld][i][j+1][k] - fsub2[fld][i][j-1][k]);
      break;
    }
    default:
    {
      printf("Error: kloop takes unknown value\n");
      break;
    }
  }
  value = value/(2.*dx);

  return value;
}

inline double deriv_yb(int fld, int i, int j, int k, int kloop_arg, FSUB_ARG)
{
  double value;

  switch(kloop_arg)
  {
    case 1:
    {
      value = (f[fld][i][INCREMENT(j)][k] - f[fld][i][DECREMENT(j)][k]);
      break;
    }
    case 2:
    case 4:
    {
      value = (fsub1[fld][i][INCREMENT(j)][k] - fsub1[fld][i][DECREMENT(j)][k]);
      break;
    }
    case 3:
    {
      value = (fsub2[fld][i][INCREMENT(j)][k] - fsub2[fld][i][DECREMENT(j)][k]);
      break;
    }
    default:
    {
      printf("Error: kloop takes unknown value\n");
      break;
    }
  }
  value = value/(2.*dx);

  return value;
}

// z-Direction
inline double deriv_z(int fld, int i, int j, int k, int kloop_arg, FSUB_ARG)
{
  double value;

  switch(kloop_arg)
  {
    case 1:
    {
      value = (f[fld][i][j][k+1] - f[fld][i][j][k-1]);
      break;
    }
    case 2:
    case 4:
    {
      value = (fsub1[fld][i][j][k+1] - fsub1[fld][i][j][k-1]);
      break;
    }
    case 3:
    {
      value = (fsub2[fld][i][j][k+1] - fsub2[fld][i][j][k-1]);
      break;
    }
    default:
    {
      printf("Error: kloop takes unknown value\n");
      break;
    }
  }
  value = value/(2.*dx);

  return value;
}

inline double deriv_zb(int fld, int i, int j, int k, int kloop_arg, FSUB_ARG)
{
  double value;

  switch(kloop_arg)
  {
    case 1:
    {
      value = (f[fld][i][j][INCREMENT(k)] - f[fld][i][j][DECREMENT(k)]);
      break;
    }
    case 2:
    case 4:
    {
      value = (fsub1[fld][i][j][INCREMENT(k)] - fsub1[fld][i][j][DECREMENT(k)]);
      break;
    }
    case 3:
    {
      value = (fsub2[fld][i][j][INCREMENT(k)] - fsub2[fld][i][j][DECREMENT(k)]);
      break;
    }
    default:
    {
      printf("Error: kloop takes unknown value\n");
      break;
    }
  }
  value = value/(2.*dx);

  return value;
}


/* Laplacian */
/*------------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------------*/
inline double lapl(int fld, int i, int j, int k, int kloop_arg, FSUB_ARG)
{
  double value;

  switch(kloop_arg)
  {
    case 1:
    {
      value = (f[fld][i][j][k+2] + f[fld][i][j][k-2]
            +f[fld][i][j+2][k] + f[fld][i][j-2][k]
            +f[fld][i+2][j][k] + f[fld][i-2][j][k]
            -6.*f[fld][i][j][k]);
      break;
    }
    case 2:
    case 4:
    {
      value = (fsub1[fld][i][j][k+2] + fsub1[fld][i][j][k-2]
            +fsub1[fld][i][j+2][k] + fsub1[fld][i][j-2][k]
            +fsub1[fld][i+2][j][k] + fsub1[fld][i-2][j][k]
            -6.*fsub1[fld][i][j][k]);
      break;
    }
    case 3:
    {
      value = (fsub2[fld][i][j][k+2] + fsub2[fld][i][j][k-2]
            +fsub2[fld][i][j+2][k] + fsub2[fld][i][j-2][k]
            +fsub2[fld][i+2][j][k] + fsub2[fld][i-2][j][k]
            -6.*fsub2[fld][i][j][k]);
      break;
    }
    default:
    {
      printf("Error: kloop takes unknown value\n");
      break;
    }
  }
  value = value/pw2(2.*dx);

  return value;
}


inline double laplb(int fld, int i, int j, int k, int kloop_arg, FSUB_ARG)
{
  double value;

  switch(kloop_arg)
  {
    case 1:
    {
      value = (f[fld][i][j][INCREMENT2(k)] + f[fld][i][j][DECREMENT2(k)]
            +f[fld][i][INCREMENT2(j)][k] + f[fld][i][DECREMENT2(j)][k]
            +f[fld][INCREMENT2(i)][j][k] + f[fld][DECREMENT2(i)][j][k]
            -6.*f[fld][i][j][k]);
      break;
    }
    case 2:
    case 4:
    {
      value = (fsub1[fld][i][j][INCREMENT2(k)] + fsub1[fld][i][j][DECREMENT2(k)]
            +fsub1[fld][i][INCREMENT2(j)][k] + fsub1[fld][i][DECREMENT2(j)][k]
            +fsub1[fld][INCREMENT2(i)][j][k] + fsub1[fld][DECREMENT2(i)][j][k]
            -6.*fsub1[fld][i][j][k]);
      break;
    }
    case 3:
    {
      value = (fsub2[fld][i][j][INCREMENT2(k)] + fsub2[fld][i][j][DECREMENT2(k)]
            +fsub2[fld][i][INCREMENT2(j)][k] + fsub2[fld][i][DECREMENT2(j)][k]
            +fsub2[fld][INCREMENT2(i)][j][k] + fsub2[fld][DECREMENT2(i)][j][k]
            -6.*fsub2[fld][i][j][k]);
      break;
    }
    default:
    {
      printf("Error: kloop takes unknown value\n");
      break;
    }
  }
  value = value/pw2(2.*dx);

  return value;
}


/* Derivation macro */
// Inflaton derivative, dir=x,y,z
#define DERIV_AX(dir) deriv_##dir(0,i,j,k,kloop_arg,FSUB)
// Inflaton derivative on boundary, dir=x,y,z
#define DERIVB_AX_Z(dir) deriv_##dir##b(0,i,j,edge_list[index],kloop_arg,FSUB)
#define DERIVB_AX_Y(dir) deriv_##dir##b(0,i,edge_list[index],j,kloop_arg,FSUB)
#define DERIVB_AX_X(dir) deriv_##dir##b(0,edge_list[index],i,j,kloop_arg,FSUB)

// Gauge derivative, dir=x,y,z, com=1~3
#define DERIV_GF(dir,com) deriv_##dir(com,i,j,k,kloop_arg,FSUB)
// Gauge derivative on boundary, dir=x,y,z, com=0~3
#define DERIVB_GF_Z(dir,com) deriv_##dir##b(com,i,j,edge_list[index],kloop_arg,FSUB)
#define DERIVB_GF_Y(dir,com) deriv_##dir##b(com,i,edge_list[index],j,kloop_arg,FSUB) 
#define DERIVB_GF_X(dir,com) deriv_##dir##b(com,edge_list[index],i,j,kloop_arg,FSUB)

// Rotation of gauge field
#define CROSSX (deriv_y(3,i,j,k,kloop_arg,FSUB) - deriv_z(2,i,j,k,kloop_arg,FSUB))
#define CROSSY (deriv_z(1,i,j,k,kloop_arg,FSUB) - deriv_x(3,i,j,k,kloop_arg,FSUB))
#define CROSSZ (deriv_x(2,i,j,k,kloop_arg,FSUB) - deriv_y(1,i,j,k,kloop_arg,FSUB))
// Rotation of gauge field on boundary, CROSSB_<direction><edge coordinate>
#define CROSSB_XZ ( deriv_yb(3,i,j,edge_list[index],kloop_arg,FSUB) - deriv_zb(2,i,j,edge_list[index],kloop_arg,FSUB) )
#define CROSSB_YZ ( deriv_zb(1,i,j,edge_list[index],kloop_arg,FSUB) - deriv_xb(3,i,j,edge_list[index],kloop_arg,FSUB) )
#define CROSSB_ZZ ( deriv_xb(2,i,j,edge_list[index],kloop_arg,FSUB) - deriv_yb(1,i,j,edge_list[index],kloop_arg,FSUB) )
#define CROSSB_XY ( deriv_yb(3,i,edge_list[index],j,kloop_arg,FSUB) - deriv_zb(2,i,edge_list[index],j,kloop_arg,FSUB) )
#define CROSSB_YY ( deriv_zb(1,i,edge_list[index],j,kloop_arg,FSUB) - deriv_xb(3,i,edge_list[index],j,kloop_arg,FSUB) )
#define CROSSB_ZY ( deriv_xb(2,i,edge_list[index],j,kloop_arg,FSUB) - deriv_yb(1,i,edge_list[index],j,kloop_arg,FSUB) )
#define CROSSB_XX ( deriv_yb(3,edge_list[index],i,j,kloop_arg,FSUB) - deriv_zb(2,edge_list[index],i,j,kloop_arg,FSUB) )
#define CROSSB_YX ( deriv_zb(1,edge_list[index],i,j,kloop_arg,FSUB) - deriv_xb(3,edge_list[index],i,j,kloop_arg,FSUB) )
#define CROSSB_ZX ( deriv_xb(2,edge_list[index],i,j,kloop_arg,FSUB) - deriv_yb(1,edge_list[index],i,j,kloop_arg,FSUB) )



/* Inflaton Gradient Energy 1/2 <|Grad(phi)|^2> = 1/2 <-phi Lapl(phi)> */
/*------------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------------*/
double gradient_energy_ax(int kloop_arg, FSUB_ARG)
{
  int i=0,j=0,k=0;
  double gradient=0.;

  switch(kloop_arg)
  {
    case 1:
    {
    #pragma omp parallel for reduction(-:gradient) num_threads(num_thread) private(j,k)
      for(i=2;i<N-2;i++)
        for(j=2;j<N-2;j++)
          for(k=2;k<N-2;k++)
            gradient -= AXION*lapl(0,i,j,k,kloop_arg,FSUB);
    #pragma omp parallel for reduction(-:gradient) num_threads(num_thread) private(j)
      for(i=0;i<N;i++)
        for(j=0;j<N;j++)
        {
          gradient -= f[0][i][j][0]*laplb(0,i,j,0,kloop_arg,FSUB) + f[0][i][j][1]*laplb(0,i,j,1,kloop_arg,FSUB) + f[0][i][j][N-2]*laplb(0,i,j,N-2,kloop_arg,FSUB) + f[0][i][j][N-1]*laplb(0,i,j,N-1,kloop_arg,FSUB); // z=0,1,N-2,N-1
          if(j<=1 || j>=N-2) continue;
          gradient -= f[0][i][0][j]*laplb(0,i,0,j,kloop_arg,FSUB) + f[0][i][1][j]*laplb(0,i,1,j,kloop_arg,FSUB) + f[0][i][N-2][j]*laplb(0,i,N-2,j,kloop_arg,FSUB) + f[0][i][N-1][j]*laplb(0,i,N-1,j,kloop_arg,FSUB); // y=0,1,N-2,N-1
          if(i<=1 || i>=N-2) continue;
          gradient -= f[0][0][i][j]*laplb(0,0,i,j,kloop_arg,FSUB) + f[0][1][i][j]*laplb(0,1,i,j,kloop_arg,FSUB) + f[0][N-2][i][j]*laplb(0,N-2,i,j,kloop_arg,FSUB) + f[0][N-1][i][j]*laplb(0,N-1,i,j,kloop_arg,FSUB); // x=0,1,N-2,N-1
        }
      gradient = gradient/pow(a,2.+2.*rescale_s);
      break;
    }
    case 2:
    case 4:
    {
    #pragma omp parallel for reduction(-:gradient) num_threads(num_thread) private(j,k)
      for(i=2;i<N-2;i++)
          for(j=2;j<N-2;j++)
            for(k=2;k<N-2;k++)
              gradient -= AXIONSUB1*lapl(0,i,j,k,kloop_arg,FSUB);
    #pragma omp parallel for reduction(-:gradient) num_threads(num_thread) private(j)
      for(i=0;i<N;i++)
        for(j=0;j<N;j++)
        {
          gradient -= fsub1[0][i][j][0]*laplb(0,i,j,0,kloop_arg,FSUB) + fsub1[0][i][j][1]*laplb(0,i,j,1,kloop_arg,FSUB) + fsub1[0][i][j][N-2]*laplb(0,i,j,N-2,kloop_arg,FSUB) + fsub1[0][i][j][N-1]*laplb(0,i,j,N-1,kloop_arg,FSUB); // z=0,1,N-2,N-1
          if(j<=1 || j>=N-2) continue;
          gradient -= fsub1[0][i][0][j]*laplb(0,i,0,j,kloop_arg,FSUB) + fsub1[0][i][1][j]*laplb(0,i,1,j,kloop_arg,FSUB) + fsub1[0][i][N-2][j]*laplb(0,i,N-2,j,kloop_arg,FSUB) + fsub1[0][i][N-1][j]*laplb(0,i,N-1,j,kloop_arg,FSUB); // y=0,1,N-2,N-1
          if(i<=1 || i>=N-2) continue;
          gradient -= fsub1[0][0][i][j]*laplb(0,0,i,j,kloop_arg,FSUB) + fsub1[0][1][i][j]*laplb(0,1,i,j,kloop_arg,FSUB) + fsub1[0][N-2][i][j]*laplb(0,N-2,i,j,kloop_arg,FSUB) + fsub1[0][N-1][i][j]*laplb(0,N-1,i,j,kloop_arg,FSUB); // x=0,1,N-2,N-1
        }
      gradient = gradient/pow(asub1,2.+2.*rescale_s);
      break;
    }
    case 3:
    {
    #pragma omp parallel for reduction(-:gradient) num_threads(num_thread) private(j,k)
      for(i=2;i<N-2;i++)
          for(j=2;j<N-2;j++)
            for(k=2;k<N-2;k++)
              gradient -= AXIONSUB2*lapl(0,i,j,k,kloop_arg,FSUB);
    #pragma omp parallel for reduction(-:gradient) num_threads(num_thread) private(j)
      for(i=0;i<N;i++)
        for(j=0;j<N;j++)
        {
          gradient -= fsub2[0][i][j][0]*laplb(0,i,j,0,kloop_arg,FSUB) + fsub2[0][i][j][1]*laplb(0,i,j,1,kloop_arg,FSUB) + fsub2[0][i][j][N-2]*laplb(0,i,j,N-2,kloop_arg,FSUB) + fsub2[0][i][j][N-1]*laplb(0,i,j,N-1,kloop_arg,FSUB); // z=0,1,N-2,N-1
          if(j<=1 || j>=N-2) continue;
          gradient -= fsub2[0][i][0][j]*laplb(0,i,0,j,kloop_arg,FSUB) + fsub2[0][i][1][j]*laplb(0,i,1,j,kloop_arg,FSUB) + fsub2[0][i][N-2][j]*laplb(0,i,N-2,j,kloop_arg,FSUB) + fsub2[0][i][N-1][j]*laplb(0,i,N-1,j,kloop_arg,FSUB); // y=0,1,N-2,N-1
          if(i<=1 || i>=N-2) continue;
          gradient -= fsub2[0][0][i][j]*laplb(0,0,i,j,kloop_arg,FSUB) + fsub2[0][1][i][j]*laplb(0,1,i,j,kloop_arg,FSUB) + fsub2[0][N-2][i][j]*laplb(0,N-2,i,j,kloop_arg,FSUB) + fsub2[0][N-1][i][j]*laplb(0,N-1,i,j,kloop_arg,FSUB); // x=0,1,N-2,N-1
        }
      gradient = gradient/pow(asub2,2.+2.*rescale_s);
      break;
    }
    default:
    {
      printf("Error: kloop takes unknown value\n");
      break;
    }
  }

  return(.5*gradient/(double)gridsize);
}


/* rho-3p (used in ad evolution) */
/*------------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------------*/
double stress_3pressure(int kloop_arg, SUB_ARG)
{
  int i,j,k;
  double rho_3p=0.;

  switch(kloop_arg)
  {
    case 1:
    {
    #pragma omp parallel for reduction(-:rho_3p) num_threads(num_thread) private(j,k)
      LOOP
      {
        rho_3p -= pw2(AXIOND);
      }
      rho_3p /= (double)gridsize*pw2(a);
      break;
    }
    case 2:
    case 4:
    {
    #pragma omp parallel for reduction(-:rho_3p) num_threads(num_thread) private(j,k)
      LOOP
      {
        rho_3p -= pw2(AXIONDSUB1);
      }
      rho_3p /= (double)gridsize*pw2(asub1);
      break;
    }
    case 3:
    {
    #pragma omp parallel for reduction(-:rho_3p) num_threads(num_thread) private(j,k)
      LOOP
      {
        rho_3p -= pw2(AXIONDSUB2);
      }
      rho_3p /= (double)gridsize*pw2(asub2);
      break;
    }
    default:
    {
      printf("Error: kloop takes unknown value\n");
      break;
    }
  }
  rho_3p += 2.*gradient_energy_ax(kloop_arg,FSUB);
  rho_3p += 4.*potential_energy(kloop_arg,FSUB);

  return rho_3p;
}

#endif