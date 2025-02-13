#ifndef _ELEMENT_OUT_
#define _ELEMENT_OUT_

#include "latticeeasy.hpp"

/* Heavily used functions */
/*------------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------------*/
int INCREMENT_OUT(int i)
{
  return( (i==N-1) ? 0 : i+1 );
}

int DECREMENT_OUT(int i)
{
  return( (i==0) ? N-1 : i-1 );
}


/* Laplacian */
/*------------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------------*/
double lapl_out(int fld, int i, int j, int k)
{
    return (f[fld][i][j][k+1] + f[fld][i][j][k-1]
            +f[fld][i][j+1][k] + f[fld][i][j-1][k]
            +f[fld][i+1][j][k] + f[fld][i-1][j][k]
            -6.*f[fld][i][j][k]);
}


double laplb_out(int fld, int i, int j, int k)
{
    return (f[fld][i][j][INCREMENT_OUT(k)] + f[fld][i][j][DECREMENT_OUT(k)]
            +f[fld][i][INCREMENT_OUT(j)][k] + f[fld][i][DECREMENT_OUT(j)][k]
            +f[fld][INCREMENT_OUT(i)][j][k] + f[fld][DECREMENT_OUT(i)][j][k]
            -6.*f[fld][i][j][k]);
}



/* Lattice Average */
/*------------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------------*/
// Field Value
/*------------------------------------------------------------------------------------*/
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
/*------------------------------------------------------------------------------------*/
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
/*------------------------------------------------------------------------------------*/
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
/*------------------------------------------------------------------------------------*/
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
/*------------------------------------------------------------------------------------*/
double field_lapl(int fld)
{
    int i,j,k;
    double value;
    double norm=1./pw2(dx);

    #pragma omp parallel for reduction(-:value) num_threads(num_thread) private(j,k)
      for(i=1;i<N-1;i++)
        for(j=1;j<N-1;j++)
          for(k=1;k<N-1;k++)
          {
            value -= FIELD(fld)*lapl_out(fld,i,j,k);
          }
    #pragma omp parallel for reduction(-:value) num_threads(num_thread) private(j)
      for(i=0;i<N;i++)
        for(j=0;j<N;j++)
        {
          value -= f[fld][i][j][0]*laplb_out(fld,i,j,0) + f[fld][i][j][N-1]*laplb_out(fld,i,j,N-1); // z=0,N-1
          if(j==0 || j==N-1) continue;
          value -= f[fld][i][0][j]*laplb_out(fld,i,0,j) + f[fld][i][N-1][j]*laplb_out(fld,i,N-1,j); // y=0,N-1
          if(i==0 || i==N-1) continue;
          value -= f[fld][0][i][j]*laplb_out(fld,0,i,j) + f[fld][N-1][i][j]*laplb_out(fld,N-1,i,j); // x=0,N-1
        }

    return value*norm/(double)gridsize;
}



// Fieldderiv*Laplacian <fd*Lapl(f)>
/*------------------------------------------------------------------------------------*/
double fieldderiv_lapl(int fld)
{
    int i,j,k;
    double value;
    double norm=1./pw2(dx);

    #pragma omp parallel for reduction(+:value) num_threads(num_thread) private(j,k)
      for(i=1;i<N-1;i++)
        for(j=1;j<N-1;j++)
          for(k=1;k<N-1;k++)
          {
            value += FIELDD(fld)*lapl_out(fld,i,j,k);
          }
    #pragma omp parallel for reduction(+:value) num_threads(num_thread) private(j)
      for(i=0;i<N;i++)
        for(j=0;j<N;j++)
        {
          value += fd[fld][i][j][0]*laplb_out(fld,i,j,0) + fd[fld][i][j][N-1]*laplb_out(fld,i,j,N-1); // z=0,N-1
          if(j==0 || j==N-1) continue;
          value += fd[fld][i][0][j]*laplb_out(fld,i,0,j) + fd[fld][i][N-1][j]*laplb_out(fld,i,N-1,j); // y=0,N-1
          if(i==0 || i==N-1) continue;
          value += fd[fld][0][i][j]*laplb_out(fld,0,i,j) + fd[fld][N-1][i][j]*laplb_out(fld,N-1,i,j); // x=0,N-1
        }

    return value*norm/(double)gridsize;
}


// Fieldderiv*dvdf
/*------------------------------------------------------------------------------------*/
double fieldderiv_dvdf(int fld)
{
    int i,j,k;
    double value;

    #pragma omp parallel for reduction(+:value) num_threads(num_thread) private(j,k)
    LOOP
        value += FIELDD(fld)*dvdf(0,i,j,k);

    return value/(double)gridsize;
}

#endif