#ifndef _ELEMENT_OUT_
#define _ELEMENT_OUT_

#include "latticeeasy.hpp"

/* Heavily used functions */
/*------------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------------*/
inline int INCREMENT_OUT(int i)
{
  return( (i==N-1) ? 0 : i+1 );
}

inline int DECREMENT_OUT(int i)
{
  return( (i==0) ? N-1 : i-1 );
}

inline int INCREMENT2_OUT(int i)
{
  return( (i>=N-2) ? ( (i==N-1) ? 1 : 0 ) : i+2 );
}

inline int DECREMENT2_OUT(int i)
{
  return( (i<=1) ? ( (i==0) ? N-2 : N-1 ) : i-2 );
}


/* Single Derivative */
/*------------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------------*/
// x-Direction
inline double deriv_out_x(int fld, int i, int j, int k)
{
  return (f[fld][i+1][j][k] - f[fld][i-1][j][k])/(2.*dx);
}

inline double deriv_out_xb(int fld, int i, int j, int k)
{
  return (f[fld][INCREMENT_OUT(i)][j][k] - f[fld][DECREMENT_OUT(i)][j][k])/(2.*dx);
}

// y-Direction
inline double deriv_out_y(int fld, int i, int j, int k)
{
  return (f[fld][i][j+1][k] - f[fld][i][j-1][k])/(2.*dx);
}

inline double deriv_out_yb(int fld, int i, int j, int k)
{
  return (f[fld][i][INCREMENT_OUT(j)][k] - f[fld][i][DECREMENT_OUT(j)][k])/(2.*dx);
}

// z-Direction
inline double deriv_out_z(int fld, int i, int j, int k)
{
  return (f[fld][i][j][k+1] - f[fld][i][j][k-1])/(2.*dx);
}

inline double deriv_out_zb(int fld, int i, int j, int k)
{
  return (f[fld][i][j][INCREMENT_OUT(k)] - f[fld][i][j][DECREMENT_OUT(k)])/(2.*dx);
}


/* Laplacian */
/*------------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------------*/
inline double lapl_out(int fld, int i, int j, int k)
{
    return (f[fld][i][j][k+2] + f[fld][i][j][k-2]
            +f[fld][i][j+2][k] + f[fld][i][j-2][k]
            +f[fld][i+2][j][k] + f[fld][i-2][j][k]
            -6.*f[fld][i][j][k])/pw2(2.*dx);
}


inline double laplb_out(int fld, int i, int j, int k)
{
    return (f[fld][i][j][INCREMENT2_OUT(k)] + f[fld][i][j][DECREMENT2_OUT(k)]
            +f[fld][i][INCREMENT2_OUT(j)][k] + f[fld][i][DECREMENT2_OUT(j)][k]
            +f[fld][INCREMENT2_OUT(i)][j][k] + f[fld][DECREMENT2_OUT(i)][j][k]
            -6.*f[fld][i][j][k])/pw2(2.*dx);
}



/* Lattice Average */
/*------------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------------*/
// Field Value
double field_value(int fld)
{
    int i,j,k;
    double value;

    #pragma omp parallel for reduction(+:value) num_threads(num_thread) private(j,k)
    LOOP
        value += FIELD(fld);

    return value/(double)gridsize;
}


// Field derivative
double field_deriv(int fld)
{
    int i,j,k;
    double value;

    #pragma omp parallel for reduction(+:value) num_threads(num_thread) private(j,k)
    LOOP
        value += FIELDD(fld);

    return value/(double)gridsize;
}


// Field Value ^2
double field_value2(int fld)
{
    int i,j,k;
    double value;

    #pragma omp parallel for reduction(+:value) num_threads(num_thread) private(j,k)
    LOOP
        value += pw2(FIELD(fld));

    return value/(double)gridsize;
}


// Field derivative ^2
double field_deriv2(int fld)
{
    int i,j,k;
    double value;

    #pragma omp parallel for reduction(+:value) num_threads(num_thread) private(j,k)
    LOOP
        value += pw2(FIELDD(fld));

    return value/(double)gridsize;
}


// Field*Laplacian <|Grad(f)|^2> = <-f Lapl(f)>
double field_lapl(int fld)
{
    int i,j,k;
    double value;

    #pragma omp parallel for reduction(-:value) num_threads(num_thread) private(j,k)
      for(i=2;i<N-2;i++)
        for(j=2;j<N-2;j++)
          for(k=2;k<N-2;k++)
          {
            value -= FIELD(fld)*lapl_out(fld,i,j,k);
          }
    #pragma omp parallel for reduction(-:value) num_threads(num_thread) private(j)
      for(i=0;i<N;i++)
        for(j=0;j<N;j++)
        {
          value -= f[fld][i][j][0]*laplb_out(fld,i,j,0) + f[fld][i][j][1]*laplb_out(fld,i,j,1) + f[fld][i][j][N-2]*laplb_out(fld,i,j,N-2) + f[fld][i][j][N-1]*laplb_out(fld,i,j,N-1); // z=0,1,N-2,N-1
          if(j<=1 || j>=N-2) continue;
          value -= f[fld][i][0][j]*laplb_out(fld,i,0,j) + f[fld][i][1][j]*laplb_out(fld,i,1,j) + f[fld][i][N-2][j]*laplb_out(fld,i,N-2,j) + f[fld][i][N-1][j]*laplb_out(fld,i,N-1,j); // y=0,1,N-2,N-1
          if(i<=1 || i>=N-2) continue;
          value -= f[fld][0][i][j]*laplb_out(fld,0,i,j) + f[fld][1][i][j]*laplb_out(fld,1,i,j) + f[fld][N-2][i][j]*laplb_out(fld,N-2,i,j) + f[fld][N-1][i][j]*laplb_out(fld,N-1,i,j); // x=0,1,N-2,N-1
        }

    return value/(double)gridsize;
}


// Fieldderiv*Laplacian <fd*Lapl(f)>
double fieldderiv_lapl(int fld)
{
    int i,j,k;
    double value;

    #pragma omp parallel for reduction(+:value) num_threads(num_thread) private(j,k)
      for(i=2;i<N-2;i++)
        for(j=2;j<N-2;j++)
          for(k=2;k<N-2;k++)
          {
            value += FIELDD(fld)*lapl_out(fld,i,j,k);
          }
    #pragma omp parallel for reduction(+:value) num_threads(num_thread) private(j)
      for(i=0;i<N;i++)
        for(j=0;j<N;j++)
        {
          value += fd[fld][i][j][0]*laplb_out(fld,i,j,0) + fd[fld][i][j][1]*laplb_out(fld,i,j,1) + fd[fld][i][j][N-2]*laplb_out(fld,i,j,N-2) + fd[fld][i][j][N-1]*laplb_out(fld,i,j,N-1); // z=0,1,N-2,N-1
          if(j<=1 || j>=N-2) continue;
          value += fd[fld][i][0][j]*laplb_out(fld,i,0,j) + fd[fld][i][1][j]*laplb_out(fld,i,1,j) + fd[fld][i][N-2][j]*laplb_out(fld,i,N-2,j) + fd[fld][i][N-1][j]*laplb_out(fld,i,N-1,j); // y=0,1,N-2,N-1
          if(i<=1 || i>=N-2) continue;
          value += fd[fld][0][i][j]*laplb_out(fld,0,i,j) + fd[fld][1][i][j]*laplb_out(fld,1,i,j) + fd[fld][N-2][i][j]*laplb_out(fld,N-2,i,j) + fd[fld][N-1][i][j]*laplb_out(fld,N-1,i,j); // x=0,1,N-2,N-1
        }

    return value/(double)gridsize;
}


// Fieldderiv*dvdf
double fieldderiv_dvdf(int fld)
{
    int i,j,k;
    double value;

    #pragma omp parallel for reduction(+:value) num_threads(num_thread) private(j,k)
    LOOP
        value += FIELDD(fld)*dvdf(i,j,k);

    return value/(double)gridsize;
}



/* Practical Function */
/*------------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------------*/
/* MACRO for practical function */
// Gauge derivative, dir=x,y,z, com=0~3
#define DERIV_OUT_GF(dir,com) deriv_out_##dir(com+1,i,j,k)
// Gauge derivative on boundary, dir=x,y,z, com=0~3
#define DERIVB_OUT_GF_Z(dir,com) deriv_out_##dir##b(com+1,i,j,edge_list[index])
#define DERIVB_OUT_GF_Y(dir,com) deriv_out_##dir##b(com+1,i,edge_list[index],j) 
#define DERIVB_OUT_GF_X(dir,com) deriv_out_##dir##b(com+1,edge_list[index],i,j)
// Rotation of gauge field
#define CROSSX_OUT (deriv_out_y(4,i,j,k) - deriv_out_z(3,i,j,k))
#define CROSSY_OUT (deriv_out_z(2,i,j,k) - deriv_out_x(4,i,j,k))
#define CROSSZ_OUT (deriv_out_x(3,i,j,k) - deriv_out_y(2,i,j,k))
// Rotation of gauge field on boundary, CROSSB_<direction><edge coordinate>
#define CROSSB_OUT_XZ ( deriv_out_yb(4,i,j,edge_list[index]) - deriv_out_zb(3,i,j,edge_list[index]) )
#define CROSSB_OUT_YZ ( deriv_out_zb(2,i,j,edge_list[index]) - deriv_out_xb(4,i,j,edge_list[index]) )
#define CROSSB_OUT_ZZ ( deriv_out_xb(3,i,j,edge_list[index]) - deriv_out_yb(2,i,j,edge_list[index]) )
#define CROSSB_OUT_XY ( deriv_out_yb(4,i,edge_list[index],j) - deriv_out_zb(3,i,edge_list[index],j) )
#define CROSSB_OUT_YY ( deriv_out_zb(2,i,edge_list[index],j) - deriv_out_xb(4,i,edge_list[index],j) )
#define CROSSB_OUT_ZY ( deriv_out_xb(3,i,edge_list[index],j) - deriv_out_yb(2,i,edge_list[index],j) )
#define CROSSB_OUT_XX ( deriv_out_yb(4,edge_list[index],i,j) - deriv_out_zb(3,edge_list[index],i,j) )
#define CROSSB_OUT_YX ( deriv_out_zb(2,edge_list[index],i,j) - deriv_out_xb(4,edge_list[index],i,j) )
#define CROSSB_OUT_ZX ( deriv_out_xb(3,edge_list[index],i,j) - deriv_out_yb(2,edge_list[index],i,j) )


// Energy density of Gauge Field
double energy_gf()
{
  int i=0,j=0,k=0;
  int edge_list[] = {0,1,N-2,N-1};
  int index;
  double energy=0.;
  double energy2=0.;
  double energy3=0.;

  // not relative boundary condition
  energy = field_deriv2(2) + field_deriv2(3) + field_deriv2(4);

  // relative boundary condition(part 1)
  #pragma omp parallel for reduction(+:energy2) num_threads(num_thread) private(j,k)
  for(i=2;i<N-2;i++)
    for(j=2;j<N-2;j++)
      for(k=2;k<N-2;k++)
        energy2 += GAUGED(1)*DERIV_OUT_GF(x,0) + GAUGED(2)*DERIV_OUT_GF(y,0) + GAUGED(3)*DERIV_OUT_GF(z,0);
  for(i=0;i<N;i++)
    for(j=0;j<N;j++)
    {
      for(index=0;index<4;index++)
      {
        // z=0,1,N-2,N-1
        energy2 += GAUGED(1)*DERIVB_OUT_GF_Z(x,0) + GAUGED(2)*DERIVB_OUT_GF_Z(y,0) + GAUGED(3)*DERIVB_OUT_GF_Z(z,0);
      }
      if(j<=1 || j>=N-2) continue;
      for(index=0;index<4;index++)
      {
        // y=0,1,N-2,N-1
        energy2 += GAUGED(1)*DERIVB_OUT_GF_Y(x,0) + GAUGED(2)*DERIVB_OUT_GF_Y(y,0) + GAUGED(3)*DERIVB_OUT_GF_Y(z,0);
      }
      if(i<=1 || i>=N-2) continue;
      for(index=0;index<4;index++)
      {
        // x=0,1,N-2,N-1
        energy2 += GAUGED(1)*DERIVB_OUT_GF_X(x,0) + GAUGED(2)*DERIVB_OUT_GF_X(y,0) + GAUGED(3)*DERIVB_OUT_GF_X(z,0);
      }
    }
  energy2 = energy2/(double)gridsize;
  energy -= energy2*2./pow(a,rescale_s);

  // relative boundary condition(part 2)
  #pragma omp parallel for reduction(+:energy3) num_threads(num_thread) private(j,k)
  for(i=2;i<N-2;i++)
    for(j=2;j<N-2;j++)
      for(k=2;k<N-2;k++)
        energy3 += pw2(DERIV_OUT_GF(x,0)) + pw2(DERIV_OUT_GF(y,0)) + pw2(DERIV_OUT_GF(z,0))
                    + pw2(DERIV_OUT_GF(y,1)) + pw2(DERIV_OUT_GF(z,2)) + pw2(DERIV_OUT_GF(x,3))
                    + pw2(DERIV_OUT_GF(z,1)) + pw2(DERIV_OUT_GF(x,2)) + pw2(DERIV_OUT_GF(y,3))
                    - 2*(DERIV_OUT_GF(y,1)*DERIV_OUT_GF(x,2) + DERIV_OUT_GF(z,2)*DERIV_OUT_GF(y,3) + DERIV_OUT_GF(x,3)*DERIV_OUT_GF(z,1));
  for(i=0;i<N;i++)
    for(j=0;j<N;j++)
    {
      for(index=0;index<4;index++)
        {
          // z=0,1,N-2,N-1
          energy3 += pw2(DERIVB_OUT_GF_Z(x,0)) + pw2(DERIVB_OUT_GF_Z(y,0)) + pw2(DERIVB_OUT_GF_Z(z,0))
                    + pw2(DERIVB_OUT_GF_Z(y,1)) + pw2(DERIVB_OUT_GF_Z(z,2)) + pw2(DERIVB_OUT_GF_Z(x,3))
                    + pw2(DERIVB_OUT_GF_Z(z,1)) + pw2(DERIVB_OUT_GF_Z(x,2)) + pw2(DERIVB_OUT_GF_Z(y,3))
                    - 2*(DERIVB_OUT_GF_Z(y,1)*DERIVB_OUT_GF_Z(x,2) + DERIVB_OUT_GF_Z(z,2)*DERIVB_OUT_GF_Z(y,3) + DERIVB_OUT_GF_Z(x,3)*DERIVB_OUT_GF_Z(z,1));
        }
        if(j<=1 || j>=N-2) continue;
        for(index=0;index<4;index++)
        {
          // y=0,1,N-2,N-1
          energy3 += pw2(DERIVB_OUT_GF_Y(x,0)) + pw2(DERIVB_OUT_GF_Y(y,0)) + pw2(DERIVB_OUT_GF_Y(z,0))
                    + pw2(DERIVB_OUT_GF_Y(y,1)) + pw2(DERIVB_OUT_GF_Y(z,2)) + pw2(DERIVB_OUT_GF_Y(x,3))
                    + pw2(DERIVB_OUT_GF_Y(z,1)) + pw2(DERIVB_OUT_GF_Y(x,2)) + pw2(DERIVB_OUT_GF_Y(y,3))
                    - 2*(DERIVB_OUT_GF_Y(y,1)*DERIVB_OUT_GF_Y(x,2) + DERIVB_OUT_GF_Y(z,2)*DERIVB_OUT_GF_Y(y,3) + DERIVB_OUT_GF_Y(x,3)*DERIVB_OUT_GF_Y(z,1));
        }
        if(i<=1 || i>=N-2) continue; // Reinterpret i as a y index. Don't double-count front and back points.
        for(index=0;index<4;index++)
        {
          // x=0,1,N-2,N-1
          energy3 += pw2(DERIVB_OUT_GF_X(x,0)) + pw2(DERIVB_OUT_GF_X(y,0)) + pw2(DERIVB_OUT_GF_X(z,0))
                    + pw2(DERIVB_OUT_GF_X(y,1)) + pw2(DERIVB_OUT_GF_X(z,2)) + pw2(DERIVB_OUT_GF_X(x,3))
                    + pw2(DERIVB_OUT_GF_X(z,1)) + pw2(DERIVB_OUT_GF_X(x,2)) + pw2(DERIVB_OUT_GF_X(y,3))
                    - 2*(DERIVB_OUT_GF_X(y,1)*DERIVB_OUT_GF_X(x,2) + DERIVB_OUT_GF_X(z,2)*DERIVB_OUT_GF_X(y,3) + DERIVB_OUT_GF_X(x,3)*DERIVB_OUT_GF_X(z,1));
        }
    }
  energy3 = energy3/(double)gridsize;
  energy += energy3/pow(a,2.*rescale_s);

  return(energy/(2.*pow(a,4.)*(double)gridsize));
}


// Slow-Roll Parameter "eta"
double eta(double gauge_energy, double slow_eps)
{
  int i,j,k;
  double deriv=0.,pot=0.,grad=0.,com=0.;
  double interaction=0.,value=0.;
  double numerator,denominator;
  int edge_list[] = {0,1,N-2,N-1};
  int index;
  double slow_eta=0.;

  deriv = field_deriv2(0);
  pot   = fieldderiv_dvdf(0);
  grad  = field_lapl(0);
  com   = fieldderiv_lapl(0);

  numerator   = -2.*deriv - 2.*pow(a,3.)*pot/ad + 4.*com/(3.*ad*pow(a,2.*rescale_s-1.)) + 2.*grad/(3.*pow(a,2.*rescale_s));
  denominator = deriv + grad/(3.*pow(a,2.*rescale_s)) + 4.*pw2(a)*gauge_energy/3.;


    for(i=2;i<N-2;i++)
        for(j=2;j<N-2;j++)
        for(k=2;k<N-2;k++)
        {
            interaction += AXIOND*((pow(a,rescale_s)*GAUGED(1)-DERIV_OUT_GF(x,0))*CROSSX_OUT + (pow(a,rescale_s)*GAUGED(2)-DERIV_OUT_GF(y,0))*CROSSY_OUT + (pow(a,rescale_s)*GAUGED(3)-DERIV_OUT_GF(z,0))*CROSSZ_OUT);
        }
    for(i=0;i<N;i++)
        for(j=0;j<N;j++)
        {
        for(index=0;index<4;index++)
        {
            value  = ( pow(a,rescale_s)*fd[2][i][j][edge_list[index]] - DERIVB_OUT_GF_Z(x,0) ) * CROSSB_OUT_XZ;
            value += ( pow(a,rescale_s)*fd[3][i][j][edge_list[index]] - DERIVB_OUT_GF_Z(y,0) ) * CROSSB_OUT_YZ;
            value += ( pow(a,rescale_s)*fd[4][i][j][edge_list[index]] - DERIVB_OUT_GF_Z(z,0) ) * CROSSB_OUT_ZZ;
            interaction += fd[0][i][j][edge_list[index]]*value;
        }
        if(j<=1 || j>=N-2) continue;
        for(index=0;index<4;index++)
        {
            value  = ( pow(a,rescale_s)*fd[2][i][edge_list[index]][j] - DERIVB_OUT_GF_Y(x,0) ) * CROSSB_OUT_XY;
            value += ( pow(a,rescale_s)*fd[3][i][edge_list[index]][j] - DERIVB_OUT_GF_Y(y,0) ) * CROSSB_OUT_YY;
            value += ( pow(a,rescale_s)*fd[4][i][edge_list[index]][j] - DERIVB_OUT_GF_Y(z,0) ) * CROSSB_OUT_ZY;
            interaction += fd[0][i][edge_list[index]][j]*value;
        }
        if(i<=1 || i>=N-2) continue;
        for(index=0;index<4;index++)
        {
            value  = ( pow(a,rescale_s)*fd[2][edge_list[index]][i][j] - DERIVB_OUT_GF_X(x,0) ) * CROSSB_OUT_XX;
            value += ( pow(a,rescale_s)*fd[3][edge_list[index]][i][j] - DERIVB_OUT_GF_X(y,0) ) * CROSSB_OUT_YX;
            value += ( pow(a,rescale_s)*fd[4][edge_list[index]][i][j] - DERIVB_OUT_GF_X(z,0) ) * CROSSB_OUT_ZX;
            interaction += fd[0][edge_list[index]][i][j]*value;
        }
        }
    interaction = interaction/(double)gridsize;

    numerator = numerator - 2.*g*interaction/(3.*ad*pow(a,2.*rescale_s+1.));
    slow_eta = numerator/denominator;
    slow_eta = slow_eta + 2.*slow_eps - 4.;

    return slow_eta;
}


// Gauge Fixing Condition
double gauge_condition()
{
  int i=0,j=0,k=0;
  int edge_list[] = {0,1,N-2,N-1};
  int index;
  double numerator=0.;
  double denominator=0.;
  double gauge = 0.;

  for(i=2;i<N-2;i++)
    for(j=2;j<N-2;j++)
      for(k=2;k<N-2;k++)
        numerator = -pow(a,rescale_s)*GAUGED(0) + DERIV_OUT_GF(x,1) + DERIV_OUT_GF(y,2) + DERIV_OUT_GF(z,3);
        denominator = sqrt( pow(a,2.*rescale_s)*pw2(GAUGED(0)) + pw2(DERIV_OUT_GF(x,1)) + pw2(DERIV_OUT_GF(y,2)) + pw2(DERIV_OUT_GF(z,3)) );
        gauge += abs(numerator/denominator);
  for(i=0;i<N;i++) // Index in x direction
    for(j=0;j<N;j++) // Index in y direction
    {
      for(index=0;index<4;index++)
        {
          // z=0,1,N-2,N-1
          numerator = -pow(a,rescale_s)*fd[1][i][j][edge_list[index]] + DERIVB_OUT_GF_Z(x,1) + DERIVB_OUT_GF_Z(y,2) + DERIVB_OUT_GF_Z(z,3);
          denominator = sqrt( pow(a,2.*rescale_s)*pw2(fd[1][i][j][edge_list[index]]) + pw2(DERIVB_OUT_GF_Z(x,1)) + pw2(DERIVB_OUT_GF_Z(y,2)) + pw2(DERIVB_OUT_GF_Z(z,3)) );
          gauge += abs(numerator/denominator);
        }
        if(j<=1 || j>=N-2) continue;
        for(index=0;index<4;index++)
        {
          // y=0,1,N-2,N-1
          numerator = -pow(a,rescale_s)*fd[1][i][edge_list[index]][j] + DERIVB_OUT_GF_Y(x,1) + DERIVB_OUT_GF_Y(y,2) + DERIVB_OUT_GF_Y(z,3);
          denominator = sqrt( pow(a,2.*rescale_s)*pw2(fd[1][i][edge_list[index]][j]) + pw2(DERIVB_OUT_GF_Y(x,1)) + pw2(DERIVB_OUT_GF_Y(y,2)) + pw2(DERIVB_OUT_GF_Y(z,3)) );
          gauge += abs(numerator/denominator);
        }
        if(i<=1 || i>=N-2) continue; // Reinterpret i as a y index. Don't double-count front and back points.
        for(index=0;index<4;index++)
        {
          // x=0,1,N-2,N-1
          numerator = -pow(a,rescale_s)*fd[1][edge_list[index]][i][j] + DERIVB_OUT_GF_X(x,1) + DERIVB_OUT_GF_X(y,2) + DERIVB_OUT_GF_X(z,3);
          denominator = sqrt( pow(a,2.*rescale_s)*pw2(fd[1][edge_list[index]][i][j]) + pw2(DERIVB_OUT_GF_X(x,1)) + pw2(DERIVB_OUT_GF_X(y,2)) + pw2(DERIVB_OUT_GF_X(z,3)) );
          gauge += abs(numerator/denominator);
        }
    }

  return (gauge/(double)gridsize);
}


// Backreaction Condition (Action for Inflaton)
double backreaction_condition()
{
  int i,j,k;
  double pot=0.;
  double interaction=0.,value=0.;
  int edge_list[] = {0,1,N-2,N-1};
  int index;

  LOOP
    pot += dvdf(i,j,k);
  pot = pot/(double)gridsize;

    #pragma omp parallel for reduction(+:interaction) num_threads(num_thread) private(j,k)
    for(i=2;i<N-2;i++)
      for(j=2;j<N-2;j++)
        for(k=2;k<N-2;k++)
        {
          interaction += ((pow(a,rescale_s)*GAUGED(1)-DERIV_OUT_GF(x,0))*CROSSX_OUT + (pow(a,rescale_s)*GAUGED(2)-DERIV_OUT_GF(y,0))*CROSSY_OUT + (pow(a,rescale_s)*GAUGED(3)-DERIV_OUT_GF(z,0))*CROSSZ_OUT);
        }
    for(i=0;i<N;i++)
      for(j=0;j<N;j++)
      {
        for(index=0;index<4;index++)
        {
          // z=0,1,N-2,N-1
          value  = ( pow(a,rescale_s)*fd[2][i][j][edge_list[index]] - DERIVB_OUT_GF_Z(x,0) ) * CROSSB_OUT_XZ;
          value += ( pow(a,rescale_s)*fd[3][i][j][edge_list[index]] - DERIVB_OUT_GF_Z(y,0) ) * CROSSB_OUT_YZ;
          value += ( pow(a,rescale_s)*fd[4][i][j][edge_list[index]] - DERIVB_OUT_GF_Z(z,0) ) * CROSSB_OUT_ZZ;
          interaction += value;
        }
        if(j<=1 || j>=N-2) continue;
        for(index=0;index<4;index++)
        {
          // y=0,1,N-2,N-1
          value  = ( pow(a,rescale_s)*fd[2][i][edge_list[index]][j] - DERIVB_OUT_GF_Y(x,0) ) * CROSSB_OUT_XY;
          value += ( pow(a,rescale_s)*fd[3][i][edge_list[index]][j] - DERIVB_OUT_GF_Y(y,0) ) * CROSSB_OUT_YY;
          value += ( pow(a,rescale_s)*fd[4][i][edge_list[index]][j] - DERIVB_OUT_GF_Y(z,0) ) * CROSSB_OUT_ZY;
          interaction += value;
        }
        if(i<=1 || i>=N-2) continue;
        for(index=0;index<4;index++)
        {
          // x=0,1,N-2,N-1
          value  = ( pow(a,rescale_s)*fd[2][edge_list[index]][i][j] - DERIVB_OUT_GF_X(x,0) ) * CROSSB_OUT_XX;
          value += ( pow(a,rescale_s)*fd[3][edge_list[index]][i][j] - DERIVB_OUT_GF_X(y,0) ) * CROSSB_OUT_YX;
          value += ( pow(a,rescale_s)*fd[4][edge_list[index]][i][j] - DERIVB_OUT_GF_X(z,0) ) * CROSSB_OUT_ZX;
          interaction += value;
        }
      }
    interaction = interaction/(double)gridsize;

    return (g*interaction/(pot*pow(a,2.*rescale_s+4.)));
}


#endif