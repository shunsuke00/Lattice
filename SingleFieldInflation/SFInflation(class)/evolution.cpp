/* This file implements time evolution function "evolve_rk4()" */

#include "latticeeasy.hpp"
#include "element_evo.hpp"


// Evolve Inflaton Field
void Field::evolve_fields(double d, int kloop_arg)
{
  DECLARE_INDICES
  int fld;

  for(fld=0;fld<nflds;fld++)
  {
    switch(kloop_arg)
    {
      case 1:
        #pragma omp parallel for num_threads(num_thread) private(j,k)
        LOOP
        {
          KFIELD(fld) = d*FIELDD(fld)/6.;
          FIELDSUB1(fld) = FIELD(fld) + d*FIELDD(fld)*0.5;
        }
        break;
      case 2:
        #pragma omp parallel for num_threads(num_thread) private(j,k)
        LOOP
        {
          KFIELD(fld) += d*FIELDDSUB1(fld)/3.;
          FIELDSUB2(fld) = FIELD(fld) + d*FIELDDSUB1(fld)*0.5;
        }
        break;
      case 3:
        #pragma omp parallel for num_threads(num_thread) private(j,k)
        LOOP
        {
          KFIELD(fld) += d*FIELDDSUB2(fld)/3.;
          FIELDSUB1(fld) = FIELD(fld) + d*FIELDDSUB2(fld);
        }
        break;
      case 4:
        #pragma omp parallel for num_threads(num_thread) private(j,k)
        LOOP
        {
          KFIELD(fld) += d*FIELDDSUB1(fld)/6.;
        }
        break;
    }
  }
}


// Evolve Time Derivative of Inflaton Field
void Field::evolve_derivs(double d, int kloop_arg)
{
  int i=0,j=0,k=0,fld;
  double laplnorm,acoeff,tempfd;

  for(fld=0;fld<nflds;fld++)
  {
    switch(kloop_arg)
    {
      case 1:
        laplnorm = 1./pw2(dx)/pow(a,2.*rescale_s);
        acoeff = ad/a;
        #pragma omp parallel for num_threads(num_thread) private(j,k,tempfd)
        for(i=1;i<N-1;i++)
          for(j=1;j<N-1;j++)
            for(k=1;k<N-1;k++)
            {
              tempfd = d*(-(rescale_s+2.)*acoeff*FIELDD(fld) + laplnorm*lapl(fld,i,j,k,kloop_arg) - pw2(a)*dvdf(fld,i,j,k,kloop_arg));
              KFIELDD(fld) = tempfd/6.;
              FIELDDSUB1(fld) = FIELDD(fld) + tempfd*0.5;
            }
        #pragma omp parallel for num_threads(num_thread) private(j,tempfd)
        for(i=0;i<N;i++)
          for(j=0;j<N;j++)
          {
            // z=0
            tempfd = d*(-(rescale_s+2.)*acoeff*fd[fld][i][j][0] + laplnorm*laplb(fld,i,j,0,kloop_arg) - pw2(a)*dvdf(fld,i,j,0,kloop_arg));
            kfd[fld][i][j][0] = tempfd/6.;
            fdsub1[fld][i][j][0] = fd[fld][i][j][0] + tempfd*0.5;
            // z=N-1
            tempfd = d*(-(rescale_s+2.)*acoeff*fd[fld][i][j][N-1] + laplnorm*laplb(fld,i,j,N-1,kloop_arg) - pw2(a)*dvdf(fld,i,j,N-1,kloop_arg));
            kfd[fld][i][j][N-1] = tempfd/6.;
            fdsub1[fld][i][j][N-1] = fd[fld][i][j][N-1] + tempfd*0.5;
            if(j==0 || j==N-1) continue;
            // y=0
            tempfd = d*(-(rescale_s+2.)*acoeff*fd[fld][i][0][j] + laplnorm*laplb(fld,i,0,j,kloop_arg) - pw2(a)*dvdf(fld,i,0,j,kloop_arg));
            kfd[fld][i][0][j] = tempfd/6.;
            fdsub1[fld][i][0][j] = fd[fld][i][0][j] + tempfd*0.5;
            // y=N-1
            tempfd = d*(-(rescale_s+2.)*acoeff*fd[fld][i][N-1][j] + laplnorm*laplb(fld,i,N-1,j,kloop_arg) - pw2(a)*dvdf(fld,i,N-1,j,kloop_arg));
            kfd[fld][i][N-1][j] = tempfd/6.;
            fdsub1[fld][i][N-1][j] = fd[fld][i][N-1][j] + tempfd*0.5;
            if(i==0 || i==N-1) continue;
            // x=0
            tempfd = d*(-(rescale_s+2.)*acoeff*fd[fld][0][i][j] + laplnorm*laplb(fld,0,i,j,kloop_arg) - pw2(a)*dvdf(fld,0,i,j,kloop_arg));
            kfd[fld][0][i][j] = tempfd/6.;
            fdsub1[fld][0][i][j] = fd[fld][0][i][j] + tempfd*0.5;
            // x=N-1
            tempfd = d*(-(rescale_s+2.)*acoeff*fd[fld][N-1][i][j] + laplnorm*laplb(fld,N-1,i,j,kloop_arg) - pw2(a)*dvdf(fld,N-1,i,j,kloop_arg));
            kfd[fld][N-1][i][j] = tempfd/6.;
            fdsub1[fld][N-1][i][j] = fd[fld][N-1][i][j] + tempfd*0.5;
          }
        break;
      case 2:
        laplnorm = 1./pw2(dx)/pow(asub1,2.*rescale_s);
        acoeff = adsub1/asub1;
        #pragma omp parallel for num_threads(num_thread) private(j,k,tempfd)
        for(i=1;i<N-1;i++)
          for(j=1;j<N-1;j++)
            for(k=1;k<N-1;k++)
            {
              tempfd = d*(-(rescale_s+2.)*acoeff*FIELDDSUB1(fld) + laplnorm*lapl(fld,i,j,k,kloop_arg) - pw2(asub1)*dvdf(fld,i,j,k,kloop_arg));
              KFIELDD(fld) += tempfd/3.;
              FIELDDSUB2(fld) = FIELDD(fld) + tempfd*0.5;
            }
        #pragma omp parallel for num_threads(num_thread) private(j,tempfd)
        for(i=0;i<N;i++)
          for(j=0;j<N;j++)
          {
            // z=0
            tempfd = d*(-(rescale_s+2.)*acoeff*fdsub1[fld][i][j][0] + laplnorm*laplb(fld,i,j,0,kloop_arg) - pw2(asub1)*dvdf(fld,i,j,0,kloop_arg));
            kfd[fld][i][j][0] += tempfd/3.;
            fdsub2[fld][i][j][0] = fd[fld][i][j][0] + tempfd*0.5;
            // z=N-1
            tempfd = d*(-(rescale_s+2.)*acoeff*fdsub1[fld][i][j][N-1] + laplnorm*laplb(fld,i,j,N-1,kloop_arg) - pw2(asub1)*dvdf(fld,i,j,N-1,kloop_arg));
            kfd[fld][i][j][N-1] += tempfd/3.;
            fdsub2[fld][i][j][N-1] = fd[fld][i][j][N-1] + tempfd*0.5;
            if(j==0 || j==N-1) continue;
            // y=0
            tempfd = d*(-(rescale_s+2.)*acoeff*fdsub1[fld][i][0][j] + laplnorm*laplb(fld,i,0,j,kloop_arg) - pw2(asub1)*dvdf(fld,i,0,j,kloop_arg));
            kfd[fld][i][0][j] += tempfd/3.;
            fdsub2[fld][i][0][j] = fd[fld][i][0][j] + tempfd*0.5;
            // y=N-1
            tempfd = d*(-(rescale_s+2.)*acoeff*fdsub1[fld][i][N-1][j] + laplnorm*laplb(fld,i,N-1,j,kloop_arg) - pw2(asub1)*dvdf(fld,i,N-1,j,kloop_arg));
            kfd[fld][i][N-1][j] += tempfd/3.;
            fdsub2[fld][i][N-1][j] = fd[fld][i][N-1][j] + tempfd*0.5;
            if(i==0 || i==N-1) continue;
            // x=0
            tempfd = d*(-(rescale_s+2.)*acoeff*fdsub1[fld][0][i][j] + laplnorm*laplb(fld,0,i,j,kloop_arg) - pw2(asub1)*dvdf(fld,0,i,j,kloop_arg));
            kfd[fld][0][i][j] += tempfd/3.;
            fdsub2[fld][0][i][j] = fd[fld][0][i][j] + tempfd*0.5;
            // x=N-1
            tempfd = d*(-(rescale_s+2.)*acoeff*fdsub1[fld][N-1][i][j] + laplnorm*laplb(fld,N-1,i,j,kloop_arg) - pw2(asub1)*dvdf(fld,N-1,i,j,kloop_arg));
            kfd[fld][N-1][i][j] += tempfd/3.;
            fdsub2[fld][N-1][i][j] = fd[fld][N-1][i][j] + tempfd*0.5;
          }
        break;
      case 3:
        laplnorm = 1./pw2(dx)/pow(asub2,2.*rescale_s);
        acoeff = adsub2/asub2;
        #pragma omp parallel for num_threads(num_thread) private(j,k,tempfd)
        for(i=1;i<N-1;i++)
          for(j=1;j<N-1;j++)
            for(k=1;k<N-1;k++)
            {
              tempfd = d*(-(rescale_s+2.)*acoeff*FIELDDSUB2(fld) + laplnorm*lapl(fld,i,j,k,kloop_arg) - pw2(asub2)*dvdf(fld,i,j,k,kloop_arg));
              KFIELDD(fld) += tempfd/3.;
              FIELDDSUB1(fld) = FIELDD(fld) + tempfd;
            }
        #pragma omp parallel for num_threads(num_thread) private(j,tempfd)
        for(i=0;i<N;i++)
          for(j=0;j<N;j++)
          {
            // z=0
            tempfd = d*(-(rescale_s+2.)*acoeff*fdsub2[fld][i][j][0] + laplnorm*laplb(fld,i,j,0,kloop_arg) - pw2(asub2)*dvdf(fld,i,j,0,kloop_arg));
            kfd[fld][i][j][0] += tempfd/3.;
            fdsub1[fld][i][j][0] = fd[fld][i][j][0] + tempfd;
            // z=N-1
            tempfd = d*(-(rescale_s+2.)*acoeff*fdsub2[fld][i][j][N-1] + laplnorm*laplb(fld,i,j,N-1,kloop_arg) - pw2(asub2)*dvdf(fld,i,j,N-1,kloop_arg));
            kfd[fld][i][j][N-1] += tempfd/3.;
            fdsub1[fld][i][j][N-1] = fd[fld][i][j][N-1] + tempfd;
            if(j==0 || j==N-1) continue;
            // y=0
            tempfd = d*(-(rescale_s+2.)*acoeff*fdsub2[fld][i][0][j] + laplnorm*laplb(fld,i,0,j,kloop_arg) - pw2(asub2)*dvdf(fld,i,0,j,kloop_arg));
            kfd[fld][i][0][j] += tempfd/3.;
            fdsub1[fld][i][0][j] = fd[fld][i][0][j] + tempfd;
            // y=N-1
            tempfd = d*(-(rescale_s+2.)*acoeff*fdsub2[fld][i][N-1][j] + laplnorm*laplb(fld,i,N-1,j,kloop_arg) - pw2(asub2)*dvdf(fld,i,N-1,j,kloop_arg));
            kfd[fld][i][N-1][j] += tempfd/3.;
            fdsub1[fld][i][N-1][j] = fd[fld][i][N-1][j] + tempfd;
            if(i==0 || i==N-1) continue;
            // x=0
            tempfd = d*(-(rescale_s+2.)*acoeff*fdsub2[fld][0][i][j] + laplnorm*laplb(fld,0,i,j,kloop_arg) - pw2(asub2)*dvdf(fld,0,i,j,kloop_arg));
            kfd[fld][0][i][j] += tempfd/3.;
            fdsub1[fld][0][i][j] = fd[fld][0][i][j] + tempfd;
            // x=N-1
            tempfd = d*(-(rescale_s+2.)*acoeff*fdsub2[fld][N-1][i][j] + laplnorm*laplb(fld,N-1,i,j,kloop_arg) - pw2(asub2)*dvdf(fld,N-1,i,j,kloop_arg));
            kfd[fld][N-1][i][j] += tempfd/3.;
            fdsub1[fld][N-1][i][j] = fd[fld][N-1][i][j] + tempfd;
          }
        break;
      case 4:
        laplnorm = 1./pw2(dx)/pow(asub1,2.*rescale_s);
        acoeff = adsub1/asub1;
        #pragma omp parallel for num_threads(num_thread) private(j,k,tempfd)
        for(i=1;i<N-1;i++)
          for(j=1;j<N-1;j++)
            for(k=1;k<N-1;k++)
            {
              KFIELDD(fld) += d*(-(rescale_s+2.)*acoeff*FIELDDSUB1(fld) + laplnorm*lapl(fld,i,j,k,kloop_arg) - pw2(asub1)*dvdf(fld,i,j,k,kloop_arg))/6.;
            }
        #pragma omp parallel for num_threads(num_thread) private(j,k,tempfd)
        for(i=0;i<N;i++)
          for(j=0;j<N;j++)
          {
            // z=0
            kfd[fld][i][j][0] += d*(-(rescale_s+2.)*acoeff*fdsub1[fld][i][j][0] + laplnorm*laplb(fld,i,j,0,kloop_arg) - pw2(asub1)*dvdf(fld,i,j,0,kloop_arg))/6.;
            // z=N-1
            kfd[fld][i][j][N-1] += d*(-(rescale_s+2.)*acoeff*fdsub1[fld][i][j][N-1] + laplnorm*laplb(fld,i,j,N-1,kloop_arg) - pw2(asub1)*dvdf(fld,i,j,N-1,kloop_arg))/6.;
            if(j==0 || j==N-1) continue;
            // y=0
            kfd[fld][i][0][j] += d*(-(rescale_s+2.)*acoeff*fdsub1[fld][i][0][j] + laplnorm*laplb(fld,i,0,j,kloop_arg) - pw2(asub1)*dvdf(fld,i,0,j,kloop_arg))/6.;
            // y=N-1
            kfd[fld][i][N-1][j] += d*(-(rescale_s+2.)*acoeff*fdsub1[fld][i][N-1][j] + laplnorm*laplb(fld,i,N-1,j,kloop_arg) - pw2(asub1)*dvdf(fld,i,N-1,j,kloop_arg))/6.;
            if(i==0 || i==N-1) continue;
            // x=0
            kfd[fld][0][i][j] += d*(-(rescale_s+2.)*acoeff*fdsub1[fld][0][i][j] + laplnorm*laplb(fld,0,i,j,kloop_arg) - pw2(asub1)*dvdf(fld,0,i,j,kloop_arg))/6.;
            // x=N-1
            kfd[fld][N-1][i][j] += d*(-(rescale_s+2.)*acoeff*fdsub1[fld][N-1][i][j] + laplnorm*laplb(fld,N-1,i,j,kloop_arg) - pw2(asub1)*dvdf(fld,N-1,i,j,kloop_arg))/6.;
          }
        break;
    }
  }
}


// Evolve Scale Factor
void Field::evolve_scale(double d, int kloop_arg)
{
  switch(kloop_arg)
  {
    case 1:
      ka = d*ad/6.;
      asub1 = a + d*ad*0.5;
      break;
    case 2:
      ka += d*adsub1/3.;
      asub2 = a + d*adsub1*0.5;
      break;
    case 3:
      ka += d*adsub2/3.;
      asub1 = a + d*adsub2;
      break;
    case 4:
      ka += d*adsub1/6.;
      break;
  }
}


// Evolve Time Derivative of Scale Factor
void Field::evolve_scaled(double d, int kloop_arg)
{
  double tempad;

  switch(kloop_arg)
  {
    case 1:
      tempad = d*(-rescale_s*pw2(ad)/a + pow(a,3.)*(stress_3pressure(kloop_arg))/(6.));
      kad = tempad/6.;
      adsub1 = ad + tempad*0.5;
      break;
    case 2:
      tempad = d*(-rescale_s*pw2(adsub1)/asub1 + pow(asub1,3.)*(stress_3pressure(kloop_arg))/(6.));
      kad += tempad/3.;
      adsub2 = ad + tempad*0.5;
      break;
    case 3:
      tempad = d*(-rescale_s*pw2(adsub2)/asub2 + pow(asub2,3.)*(stress_3pressure(kloop_arg))/(6.));
      kad += tempad/3.;
      adsub1 = ad + tempad;
      break;
    case 4:
      kad += d*(-rescale_s*pw2(adsub1)/asub1 + pow(asub1,3.)*(stress_3pressure(kloop_arg))/(6.))/6.;
      break;
  }
}


// Important Function of Evolve 1 Step
void Field::evolve_rk4(double dt)
{
  DECLARE_INDICES
  int kloop,fld;

  /* calculate delta */
  for(kloop=1;kloop<=4;kloop++)
  {
    evolve_fields(dt,kloop);
    evolve_derivs(dt,kloop);
    evolve_scale (dt,kloop);
    evolve_scaled(dt,kloop);
  }

  /* evolution */
  for(fld=0;fld<nflds;fld++)
  {
    LOOP
    {
      FIELD(fld)  += KFIELD(fld);
      FIELDD(fld) += KFIELDD(fld);
    }
  }
  a  += ka;
  ad += kad;

  t += dt;
}
