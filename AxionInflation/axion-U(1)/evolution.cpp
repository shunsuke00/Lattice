/* This file implements time evolution function "evolve_rk4()" */

#include "latticeeasy.hpp"
#include "element_evo.hpp"


// Evolve Inflaton(Axion) and Gauge Field
void evolve_fields(double d, int kloop_arg, SUB_ARG, double kf[nflds][N][N][N])
{
  DECLARE_INDICES
  int fld;

  for(fld=0;fld<nflds;fld++)
  {
    switch(kloop_arg)
    {
      case 1:
      {
        #pragma omp parallel for num_threads(num_thread) private(j,k)
        LOOP
        {
          KFIELD(fld) = d*FIELDD(fld)/6.;
          FIELDSUB1(fld) = FIELD(fld) + d*FIELDD(fld)/2.;
        }
        break;
      }
      case 2:
      {
        #pragma omp parallel for num_threads(num_thread) private(j,k)
        LOOP
        {
          KFIELD(fld) += d*FIELDDSUB1(fld)/3.;
          FIELDSUB2(fld) = FIELD(fld) + d*FIELDDSUB1(fld)/2.;
        }
        break;
      }
      case 3:
      {
        #pragma omp parallel for num_threads(num_thread) private(j,k)
        LOOP
        {
          KFIELD(fld) += d*FIELDDSUB2(fld)/3.;
          FIELDSUB1(fld) = FIELD(fld) + d*FIELDDSUB2(fld);
        }
        break;
      }
      case 4:
      {
        #pragma omp parallel for num_threads(num_thread) private(j,k)
        LOOP
        {
          KFIELD(fld) += d*FIELDDSUB1(fld)/6.;
        }
        break;
      }
      default:
      {
        printf("Error: kloop takes unknown value\n");
        break;
      }
    }
  }
}


// Evolve Time Derivative of Inflaton(Axion) Field
void evolve_derivs_ax(double d, int kloop_arg, SUB_ARG, double kfd[nflds][N][N][N])
{
  int i=0,j=0,k=0;
  double laplnorm,backnorm,temp,back_reaction;
  int edge_list[] = {0,1,N-2,N-1};
  int index;

  switch(kloop_arg)
  {
    case 1:
    {
      laplnorm = 1./pow(a,2.*rescale_s);
      backnorm = 1./pow(a,2.+2.*rescale_s);
      #pragma omp parallel for num_threads(num_thread) private(j,k,temp,back_reaction)
      for(i=2;i<N-2;i++)
        for(j=2;j<N-2;j++)
          for(k=2;k<N-2;k++)
          {
            temp = d*(-(rescale_s+2.)*ad*AXIOND/a + laplnorm*lapl(0,i,j,k,kloop_arg,FSUB) - pw2(a)*dvdf(i,j,k,kloop_arg,FSUB));
            back_reaction  = ( pow(a,rescale_s)*GAUGED(1) - DERIV_GF(x,0) ) * CROSSX;
            back_reaction += ( pow(a,rescale_s)*GAUGED(2) - DERIV_GF(y,0) ) * CROSSY;
            back_reaction += ( pow(a,rescale_s)*GAUGED(3) - DERIV_GF(z,0) ) * CROSSZ;
            temp -= d*g*backnorm*back_reaction;
            KAXIOND = temp/6.;
            AXIONDSUB1 = AXIOND + temp/2.;
          }
      #pragma omp parallel for num_threads(num_thread) private(j,temp,back_reaction,index)
      for(i=0;i<N;i++) // Index in x direction
        for(j=0;j<N;j++) // Index in y direction
        {
          for(index=0;index<4;index++)
          {
            // z=0,1,N-2,N-1
            temp = d*(-(rescale_s+2.)*ad*fd[0][i][j][edge_list[index]]/a + laplnorm*laplb(0,i,j,edge_list[index],kloop_arg,FSUB) - pw2(a)*dvdf(i,j,edge_list[index],kloop_arg,FSUB));
            back_reaction  = ( pow(a,rescale_s)*fd[2][i][j][edge_list[index]] - DERIVB_GF_Z(x,0) ) * CROSSB_XZ;
            back_reaction += ( pow(a,rescale_s)*fd[3][i][j][edge_list[index]] - DERIVB_GF_Z(y,0) ) * CROSSB_YZ;
            back_reaction += ( pow(a,rescale_s)*fd[4][i][j][edge_list[index]] - DERIVB_GF_Z(z,0) ) * CROSSB_ZZ;
            temp -= d*g*backnorm*back_reaction;
            kfd[0][i][j][edge_list[index]] = temp/6.;
            fdsub1[0][i][j][edge_list[index]] = fd[0][i][j][edge_list[index]] + temp/2.;
          }
          if(j<=1 || j>=N-2) continue;
          for(index=0;index<4;index++)
          {
            // y=0,1,N-2,N-1
            temp = d*(-(rescale_s+2.)*ad*fd[0][i][edge_list[index]][j]/a + laplnorm*laplb(0,i,edge_list[index],j,kloop_arg,FSUB) - pw2(a)*dvdf(i,edge_list[index],j,kloop_arg,FSUB));
            back_reaction  = ( pow(a,rescale_s)*fd[2][i][edge_list[index]][j] - DERIVB_GF_Y(x,0) ) * CROSSB_XY;
            back_reaction += ( pow(a,rescale_s)*fd[3][i][edge_list[index]][j] - DERIVB_GF_Y(y,0) ) * CROSSB_YY;
            back_reaction += ( pow(a,rescale_s)*fd[4][i][edge_list[index]][j] - DERIVB_GF_Y(z,0) ) * CROSSB_ZY;
            temp -= d*g*backnorm*back_reaction;
            kfd[0][i][edge_list[index]][j] = temp/6.;
            fdsub1[0][i][edge_list[index]][j] = fd[0][i][edge_list[index]][j] + temp/2.;
          }
          if(i<=1 || i>=N-2) continue; // Reinterpret i as a y index. Don't double-count front and back points.
          for(index=0;index<4;index++)
          {
            // x=0,1,N-2,N-1
            temp = d*(-(rescale_s+2.)*ad*fd[0][edge_list[index]][i][j]/a + laplnorm*laplb(0,edge_list[index],i,j,kloop_arg,FSUB) - pw2(a)*dvdf(edge_list[index],i,j,kloop_arg,FSUB));
            back_reaction  = ( pow(a,rescale_s)*fd[2][edge_list[index]][i][j] - DERIVB_GF_X(x,0) ) * CROSSB_XX;
            back_reaction += ( pow(a,rescale_s)*fd[3][edge_list[index]][i][j] - DERIVB_GF_X(y,0) ) * CROSSB_YX;
            back_reaction += ( pow(a,rescale_s)*fd[4][edge_list[index]][i][j] - DERIVB_GF_X(z,0) ) * CROSSB_ZX;
            temp -= d*g*backnorm*back_reaction;
            kfd[0][edge_list[index]][i][j] = temp/6.;
            fdsub1[0][edge_list[index]][i][j] = fd[0][edge_list[index]][i][j] + temp/2.;
          }
        }
      break;
    }
    case 2:
    {
      laplnorm = 1./pow(asub1,2.*rescale_s);
      backnorm = 1./pow(asub1,2.+2.*rescale_s);
      #pragma omp parallel for num_threads(num_thread) private(j,k,temp,back_reaction)
      for(i=2;i<N-2;i++)
        for(j=2;j<N-2;j++)
          for(k=2;k<N-2;k++)
          {
            temp = d*(-(rescale_s+2.)*adsub1*AXIONDSUB1/asub1 + laplnorm*lapl(0,i,j,k,kloop_arg,FSUB) - pw2(asub1)*dvdf(i,j,k,kloop_arg,FSUB));
            back_reaction  = ( pow(asub1,rescale_s)*GAUGEDSUB1(1) - DERIV_GF(x,0) ) * CROSSX;
            back_reaction += ( pow(asub1,rescale_s)*GAUGEDSUB1(2) - DERIV_GF(y,0) ) * CROSSY;
            back_reaction += ( pow(asub1,rescale_s)*GAUGEDSUB1(3) - DERIV_GF(z,0) ) * CROSSZ;
            temp -= d*g*backnorm*back_reaction;
            KAXIOND += temp/3.;
            AXIONDSUB2 = AXIOND + temp/2.;
          }
      #pragma omp parallel for num_threads(num_thread) private(j,temp,back_reaction,index)
      for(i=0;i<N;i++) // Index in x direction
        for(j=0;j<N;j++) // Index in y direction
        {
          for(index=0;index<4;index++)
          {
            // z=0,1,N-2,N-1
            temp = d*(-(rescale_s+2.)*adsub1*fdsub1[0][i][j][edge_list[index]]/asub1 + laplnorm*laplb(0,i,j,edge_list[index],kloop_arg,FSUB) - pw2(asub1)*dvdf(i,j,edge_list[index],kloop_arg,FSUB));
            back_reaction  = ( pow(asub1,rescale_s)*fdsub1[2][i][j][edge_list[index]] - DERIVB_GF_Z(x,0) ) * CROSSB_XZ;
            back_reaction += ( pow(asub1,rescale_s)*fdsub1[3][i][j][edge_list[index]] - DERIVB_GF_Z(y,0) ) * CROSSB_YZ;
            back_reaction += ( pow(asub1,rescale_s)*fdsub1[4][i][j][edge_list[index]] - DERIVB_GF_Z(z,0) ) * CROSSB_ZZ;
            temp -= d*g*backnorm*back_reaction;
            kfd[0][i][j][edge_list[index]] += temp/3.;
            fdsub2[0][i][j][edge_list[index]] = fd[0][i][j][edge_list[index]] + temp/2.;
          }
          if(j<=1 || j>=N-2) continue;
          for(index=0;index<4;index++)
          {
            // y=0,1,N-2,N-1
            temp = d*(-(rescale_s+2.)*adsub1*fdsub1[0][i][edge_list[index]][j]/asub1 + laplnorm*laplb(0,i,edge_list[index],j,kloop_arg,FSUB) - pw2(asub1)*dvdf(i,edge_list[index],j,kloop_arg,FSUB));
            back_reaction  = ( pow(asub1,rescale_s)*fdsub1[2][i][edge_list[index]][j] - DERIVB_GF_Y(x,0) ) * CROSSB_XY;
            back_reaction += ( pow(asub1,rescale_s)*fdsub1[3][i][edge_list[index]][j] - DERIVB_GF_Y(y,0) ) * CROSSB_YY;
            back_reaction += ( pow(asub1,rescale_s)*fdsub1[4][i][edge_list[index]][j] - DERIVB_GF_Y(z,0) ) * CROSSB_ZY;
            temp -= d*g*backnorm*back_reaction;
            kfd[0][i][edge_list[index]][j] += temp/3.;
            fdsub2[0][i][edge_list[index]][j] = fd[0][i][edge_list[index]][j] + temp/2.;
          }
          if(i<=1 || i>=N-2) continue; // Reinterpret i as a y index. Don't double-count front and back points.
          for(index=0;index<4;index++)
          {
            // x=0,1,N-2,N-1
            temp = d*(-(rescale_s+2.)*adsub1*fdsub1[0][edge_list[index]][i][j]/asub1 + laplnorm*laplb(0,edge_list[index],i,j,kloop_arg,FSUB) - pw2(asub1)*dvdf(edge_list[index],i,j,kloop_arg,FSUB));
            back_reaction  = ( pow(asub1,rescale_s)*fdsub1[2][edge_list[index]][i][j] - DERIVB_GF_X(x,0) ) * CROSSB_XX;
            back_reaction += ( pow(asub1,rescale_s)*fdsub1[3][edge_list[index]][i][j] - DERIVB_GF_X(y,0) ) * CROSSB_YX;
            back_reaction += ( pow(asub1,rescale_s)*fdsub1[4][edge_list[index]][i][j] - DERIVB_GF_X(z,0) ) * CROSSB_ZX;
            temp -= d*g*backnorm*back_reaction;
            kfd[0][edge_list[index]][i][j] += temp/3.;
            fdsub2[0][edge_list[index]][i][j] = fd[0][edge_list[index]][i][j] + temp/2.;
          }
        }
      break;
    }
    case 3:
    {
      laplnorm = 1./pow(asub2,2.*rescale_s);
      backnorm = 1./pow(asub2,2.+2.*rescale_s);
      #pragma omp parallel for num_threads(num_thread) private(j,k,temp,back_reaction)
      for(i=2;i<N-2;i++)
        for(j=2;j<N-2;j++)
          for(k=2;k<N-2;k++)
          {
            temp = d*(-(rescale_s+2.)*adsub2*AXIONDSUB2/asub2 + laplnorm*lapl(0,i,j,k,kloop_arg,FSUB) - pw2(asub2)*dvdf(i,j,k,kloop_arg,FSUB));
            back_reaction  = ( pow(asub2,rescale_s)*GAUGEDSUB2(1) - DERIV_GF(x,0) ) * CROSSX;
            back_reaction += ( pow(asub2,rescale_s)*GAUGEDSUB2(2) - DERIV_GF(y,0) ) * CROSSY;
            back_reaction += ( pow(asub2,rescale_s)*GAUGEDSUB2(3) - DERIV_GF(z,0) ) * CROSSZ;
            temp -= d*g*backnorm*back_reaction;
            KAXIOND += temp/3.;
            AXIONDSUB1 = AXIOND + temp;
          }
      #pragma omp parallel for num_threads(num_thread) private(j,temp,back_reaction,index)
      for(i=0;i<N;i++) // Index in x direction
        for(j=0;j<N;j++) // Index in y direction
        {
          for(index=0;index<4;index++)
          {
            // z=0,1,N-2,N-1
            temp = d*(-(rescale_s+2.)*adsub2*fdsub2[0][i][j][edge_list[index]]/asub2 + laplnorm*laplb(0,i,j,edge_list[index],kloop_arg,FSUB) - pw2(asub2)*dvdf(i,j,edge_list[index],kloop_arg,FSUB));
            back_reaction  = ( pow(asub2,rescale_s)*fdsub2[2][i][j][edge_list[index]] - DERIVB_GF_Z(x,0) ) * CROSSB_XZ;
            back_reaction += ( pow(asub2,rescale_s)*fdsub2[3][i][j][edge_list[index]] - DERIVB_GF_Z(y,0) ) * CROSSB_YZ;
            back_reaction += ( pow(asub2,rescale_s)*fdsub2[4][i][j][edge_list[index]] - DERIVB_GF_Z(z,0) ) * CROSSB_ZZ;
            temp -= d*g*backnorm*back_reaction;
            kfd[0][i][j][edge_list[index]] += temp/3.;
            fdsub1[0][i][j][edge_list[index]] = fd[0][i][j][edge_list[index]] + temp;
          }
          if(j<=1 || j>=N-2) continue;
          for(index=0;index<4;index++)
          {
            // y=0,1,N-2,N-1
            temp = d*(-(rescale_s+2.)*adsub2*fdsub2[0][i][edge_list[index]][j]/asub2 + laplnorm*laplb(0,i,edge_list[index],j,kloop_arg,FSUB) - pw2(asub2)*dvdf(i,edge_list[index],j,kloop_arg,FSUB));
            back_reaction  = ( pow(asub2,rescale_s)*fdsub2[2][i][edge_list[index]][j] - DERIVB_GF_Y(x,0) ) * CROSSB_XY;
            back_reaction += ( pow(asub2,rescale_s)*fdsub2[3][i][edge_list[index]][j] - DERIVB_GF_Y(y,0) ) * CROSSB_YY;
            back_reaction += ( pow(asub2,rescale_s)*fdsub2[4][i][edge_list[index]][j] - DERIVB_GF_Y(z,0) ) * CROSSB_ZY;
            temp -= d*g*backnorm*back_reaction;
            kfd[0][i][edge_list[index]][j] += temp/3.;
            fdsub1[0][i][edge_list[index]][j] = fd[0][i][edge_list[index]][j] + temp;
          }
          if(i<=1 || i>=N-2) continue; // Reinterpret i as a y index. Don't double-count front and back points.
          for(index=0;index<4;index++)
          {
            // x=0,1,N-2,N-1
            temp = d*(-(rescale_s+2.)*adsub2*fdsub2[0][edge_list[index]][i][j]/asub2 + laplnorm*laplb(0,edge_list[index],i,j,kloop_arg,FSUB) - pw2(asub2)*dvdf(edge_list[index],i,j,kloop_arg,FSUB));
            back_reaction  = ( pow(asub2,rescale_s)*fdsub2[2][edge_list[index]][i][j] - DERIVB_GF_X(x,0) ) * CROSSB_XX;
            back_reaction += ( pow(asub2,rescale_s)*fdsub2[3][edge_list[index]][i][j] - DERIVB_GF_X(y,0) ) * CROSSB_YX;
            back_reaction += ( pow(asub2,rescale_s)*fdsub2[4][edge_list[index]][i][j] - DERIVB_GF_X(z,0) ) * CROSSB_ZX;
            temp -= d*g*backnorm*back_reaction;
            kfd[0][edge_list[index]][i][j] += temp/3.;
            fdsub1[0][edge_list[index]][i][j] = fd[0][edge_list[index]][i][j] + temp;
          }
        }
      break;
    }
    case 4:
    {
      laplnorm = 1./pow(asub1,2.*rescale_s);
      backnorm = 1./pow(asub1,2.+2.*rescale_s);
      #pragma omp parallel for num_threads(num_thread) private(j,k,temp,back_reaction)
      for(i=2;i<N-2;i++)
        for(j=2;j<N-2;j++)
          for(k=2;k<N-2;k++)
          {
            temp = d*(-(rescale_s+2.)*adsub1*AXIONDSUB1/asub1 + laplnorm*lapl(0,i,j,k,kloop_arg,FSUB) - pw2(asub1)*dvdf(i,j,k,kloop_arg,FSUB));
            back_reaction  = ( pow(asub1,rescale_s)*GAUGEDSUB1(1) - DERIV_GF(x,0) ) * CROSSX;
            back_reaction += ( pow(asub1,rescale_s)*GAUGEDSUB1(2) - DERIV_GF(y,0) ) * CROSSY;
            back_reaction += ( pow(asub1,rescale_s)*GAUGEDSUB1(3) - DERIV_GF(z,0) ) * CROSSZ;
            temp -= d*g*backnorm*back_reaction;
            KAXIOND += temp/6.;
          }
      #pragma omp parallel for num_threads(num_thread) private(j,temp,back_reaction,index)
      for(i=0;i<N;i++) // Index in x direction
        for(j=0;j<N;j++) // Index in y direction
        {
          for(index=0;index<4;index++)
          {
            // z=0,1,N-2,N-1
            temp = d*(-(rescale_s+2.)*adsub1*fdsub1[0][i][j][edge_list[index]]/asub1 + laplnorm*laplb(0,i,j,edge_list[index],kloop_arg,FSUB) - pw2(asub1)*dvdf(i,j,edge_list[index],kloop_arg,FSUB));
            back_reaction  = ( pow(asub1,rescale_s)*fdsub1[2][i][j][edge_list[index]] - DERIVB_GF_Z(x,0) ) * CROSSB_XZ;
            back_reaction += ( pow(asub1,rescale_s)*fdsub1[3][i][j][edge_list[index]] - DERIVB_GF_Z(y,0) ) * CROSSB_YZ;
            back_reaction += ( pow(asub1,rescale_s)*fdsub1[4][i][j][edge_list[index]] - DERIVB_GF_Z(z,0) ) * CROSSB_ZZ;
            temp -= d*g*backnorm*back_reaction;
            kfd[0][i][j][edge_list[index]] += temp/6.;
          }
          if(j<=1 || j>=N-2) continue;
          for(index=0;index<4;index++)
          {
            // y=0,1,N-2,N-1
            temp = d*(-(rescale_s+2.)*adsub1*fdsub1[0][i][edge_list[index]][j]/asub1 + laplnorm*laplb(0,i,edge_list[index],j,kloop_arg,FSUB) - pw2(asub1)*dvdf(i,edge_list[index],j,kloop_arg,FSUB));
            back_reaction  = ( pow(asub1,rescale_s)*fdsub1[2][i][edge_list[index]][j] - DERIVB_GF_Y(x,0) ) * CROSSB_XY;
            back_reaction += ( pow(asub1,rescale_s)*fdsub1[3][i][edge_list[index]][j] - DERIVB_GF_Y(y,0) ) * CROSSB_YY;
            back_reaction += ( pow(asub1,rescale_s)*fdsub1[4][i][edge_list[index]][j] - DERIVB_GF_Y(z,0) ) * CROSSB_ZY;
            temp -= d*g*backnorm*back_reaction;
            kfd[0][i][edge_list[index]][j] += temp/6.;
          }
          if(i<=1 || i>=N-2) continue; // Reinterpret i as a y index. Don't double-count front and back points.
          for(index=0;index<4;index++)
          {
            // x=0,1,N-2,N-1
            temp = d*(-(rescale_s+2.)*adsub1*fdsub1[0][edge_list[index]][i][j]/asub1 + laplnorm*laplb(0,edge_list[index],i,j,kloop_arg,FSUB) - pw2(asub1)*dvdf(edge_list[index],i,j,kloop_arg,FSUB));
            back_reaction  = ( pow(asub1,rescale_s)*fdsub1[2][edge_list[index]][i][j] - DERIVB_GF_X(x,0) ) * CROSSB_XX;
            back_reaction += ( pow(asub1,rescale_s)*fdsub1[3][edge_list[index]][i][j] - DERIVB_GF_X(y,0) ) * CROSSB_YX;
            back_reaction += ( pow(asub1,rescale_s)*fdsub1[4][edge_list[index]][i][j] - DERIVB_GF_X(z,0) ) * CROSSB_ZX;
            temp -= d*g*backnorm*back_reaction;
            kfd[0][edge_list[index]][i][j] += temp/6.;
          }
        }
      break;
    }
    default:
    {
      printf("Error: kloop takes unknown value\n");
      break;
    }
  }
}


// Evolve Time Derivative of Gauge Field (0 component)
void evolve_derivs_gf0(double d, int kloop_arg, SUB_ARG, double kfd[nflds][N][N][N])
{
  int i=0,j=0,k=0;
  double laplnorm,temp,interaction;
  int edge_list[] = {0,1,N-2,N-1};
  int index;

  switch(kloop_arg)
  {
    case 1:
    {
      laplnorm = 1./pow(a,2.*rescale_s);
      #pragma omp parallel for num_threads(num_thread) private(j,k,temp,interaction)
      for(i=2;i<N-2;i++)
        for(j=2;j<N-2;j++)
          for(k=2;k<N-2;k++)
          {
            temp = d*(-rescale_s*ad*GAUGED(0)/a + laplnorm*lapl(1,i,j,k,kloop_arg,FSUB));
            interaction  = DERIV_AX(x)*CROSSX + DERIV_AX(y)*CROSSY + DERIV_AX(z)*CROSSZ;
            temp += d*g*laplnorm*interaction;
            KGAUGED(0) = temp/6.;
            GAUGEDSUB1(0) = GAUGED(0) + temp/2.;
          }
      #pragma omp parallel for num_threads(num_thread) private(j,temp,interaction,index)
      for(i=0;i<N;i++) // Index in x direction
        for(j=0;j<N;j++) // Index in y direction
        {
          for(index=0;index<4;index++)
          {
            // z=0,1,N-2,N-1
            temp = d*(-rescale_s*ad*fd[1][i][j][edge_list[index]]/a + laplnorm*laplb(1,i,j,edge_list[index],kloop_arg,FSUB));
            interaction  = DERIVB_AX_Z(x)*CROSSB_XZ + DERIVB_AX_Z(y)*CROSSB_YZ + DERIVB_AX_Z(z)*CROSSB_ZZ;
            temp += d*g*laplnorm*interaction;
            kfd[1][i][j][edge_list[index]] = temp/6.;
            fdsub1[1][i][j][edge_list[index]] = fd[1][i][j][edge_list[index]] + temp/2.;
          }
          if(j<=1 || j>=N-2) continue;
          for(index=0;index<4;index++)
          {
            // y=0,1,N-2,N-1
            temp = d*(-rescale_s*ad*fd[1][i][edge_list[index]][j]/a + laplnorm*laplb(1,i,edge_list[index],j,kloop_arg,FSUB));
            interaction  = DERIVB_AX_Y(x)*CROSSB_XY + DERIVB_AX_Y(y)*CROSSB_YY + DERIVB_AX_Y(z)*CROSSB_ZY;
            temp += d*g*laplnorm*interaction;
            kfd[1][i][edge_list[index]][j] = temp/6.;
            fdsub1[1][i][edge_list[index]][j] = fd[1][i][edge_list[index]][j] + temp/2.;
          }
          if(i<=1 || i>=N-2) continue;
          for(index=0;index<4;index++)
          {
            // x=0,1,N-2,N-1
            temp = d*(-rescale_s*ad*fd[1][edge_list[index]][i][j]/a + laplnorm*laplb(1,edge_list[index],i,j,kloop_arg,FSUB));
            interaction  = DERIVB_AX_X(x)*CROSSB_XX + DERIVB_AX_X(y)*CROSSB_YX + DERIVB_AX_X(z)*CROSSB_ZX;
            temp += d*g*laplnorm*interaction;
            kfd[1][edge_list[index]][i][j] = temp/6.;
            fdsub1[1][edge_list[index]][i][j] = fd[1][edge_list[index]][i][j] + temp/2.;
          }
        }
      break;
    }
    case 2:
    {
      laplnorm = 1./pow(asub1,2.*rescale_s);
      #pragma omp parallel for num_threads(num_thread) private(j,k,temp,interaction)
      for(i=2;i<N-2;i++)
        for(j=2;j<N-2;j++)
          for(k=2;k<N-2;k++)
          {
            temp = d*(-rescale_s*adsub1*GAUGEDSUB1(0)/asub1 + laplnorm*lapl(1,i,j,k,kloop_arg,FSUB));
            interaction  = DERIV_AX(x)*CROSSX + DERIV_AX(y)*CROSSY + DERIV_AX(z)*CROSSZ;
            temp += d*g*laplnorm*interaction;
            KGAUGED(0) += temp/3.;
            GAUGEDSUB2(0) = GAUGED(0) + temp/2.;
          }
      #pragma omp parallel for num_threads(num_thread) private(j,temp,interaction,index)
      for(i=0;i<N;i++) // Index in x direction
        for(j=0;j<N;j++) // Index in y direction
        {
          for(index=0;index<4;index++)
          {
            // z=0,1,N-2,N-1
            temp = d*(-rescale_s*adsub1*fdsub1[1][i][j][edge_list[index]]/asub1 + laplnorm*laplb(1,i,j,edge_list[index],kloop_arg,FSUB));
            interaction  = DERIVB_AX_Z(x)*CROSSB_XZ + DERIVB_AX_Z(y)*CROSSB_YZ + DERIVB_AX_Z(z)*CROSSB_ZZ;
            temp += d*g*laplnorm*interaction;
            kfd[1][i][j][edge_list[index]] += temp/3.;
            fdsub2[1][i][j][edge_list[index]] = fd[1][i][j][edge_list[index]] + temp/2.;
          }
          if(j<=1 || j>=N-2) continue;
          for(index=0;index<4;index++)
          {
            // y=0,1,N-2,N-1
            temp = d*(-rescale_s*adsub1*fdsub1[1][i][edge_list[index]][j]/asub1 + laplnorm*laplb(1,i,edge_list[index],j,kloop_arg,FSUB));
            interaction  = DERIVB_AX_Y(x)*CROSSB_XY + DERIVB_AX_Y(y)*CROSSB_YY + DERIVB_AX_Y(z)*CROSSB_ZY;
            temp += d*g*laplnorm*interaction;
            kfd[1][i][edge_list[index]][j] += temp/3.;
            fdsub2[1][i][edge_list[index]][j] = fd[1][i][edge_list[index]][j] + temp/2.;
          }
          if(i<=1 || i>=N-2) continue;
          for(index=0;index<4;index++)
          {
            // x=0,1,N-2,N-1
            temp = d*(-rescale_s*adsub1*fdsub1[1][edge_list[index]][i][j]/asub1 + laplnorm*laplb(1,edge_list[index],i,j,kloop_arg,FSUB));
            interaction  = DERIVB_AX_X(x)*CROSSB_XX + DERIVB_AX_X(y)*CROSSB_YX + DERIVB_AX_X(z)*CROSSB_ZX;
            temp += d*g*laplnorm*interaction;
            kfd[1][edge_list[index]][i][j] += temp/3.;
            fdsub2[1][edge_list[index]][i][j] = fd[1][edge_list[index]][i][j] + temp/2.;
          }
        }
      break;
    }
    case 3:
    {
      laplnorm = 1./pow(asub2,2.*rescale_s);
      #pragma omp parallel for num_threads(num_thread) private(j,k,temp,interaction)
      for(i=2;i<N-2;i++)
        for(j=2;j<N-2;j++)
          for(k=2;k<N-2;k++)
          {
            temp = d*(-rescale_s*adsub2*GAUGEDSUB2(0)/asub2 + laplnorm*lapl(1,i,j,k,kloop_arg,FSUB));
            interaction  = DERIV_AX(x)*CROSSX + DERIV_AX(y)*CROSSY + DERIV_AX(z)*CROSSZ;
            temp += d*g*laplnorm*interaction;
            KGAUGED(0) += temp/3.;
            GAUGEDSUB1(0) = GAUGED(0) + temp;
          }
      #pragma omp parallel for num_threads(num_thread) private(j,temp,interaction,index)
      for(i=0;i<N;i++) // Index in x direction
        for(j=0;j<N;j++) // Index in y direction
        {
          for(index=0;index<4;index++)
          {
            // z=0,1,N-2,N-1
            temp = d*(-rescale_s*adsub2*fdsub2[1][i][j][edge_list[index]]/asub2 + laplnorm*laplb(1,i,j,edge_list[index],kloop_arg,FSUB));
            interaction  = DERIVB_AX_Z(x)*CROSSB_XZ + DERIVB_AX_Z(y)*CROSSB_YZ + DERIVB_AX_Z(z)*CROSSB_ZZ;
            temp += d*g*laplnorm*interaction;
            kfd[1][i][j][edge_list[index]] += temp/3.;
            fdsub1[1][i][j][edge_list[index]] = fd[1][i][j][edge_list[index]] + temp;
          }
          if(j<=1 || j>=N-2) continue;
          for(index=0;index<4;index++)
          {
            // y=0,1,N-2,N-1
            temp = d*(-rescale_s*adsub2*fdsub2[1][i][edge_list[index]][j]/asub2 + laplnorm*laplb(1,i,edge_list[index],j,kloop_arg,FSUB));
            interaction  = DERIVB_AX_Y(x)*CROSSB_XY + DERIVB_AX_Y(y)*CROSSB_YY + DERIVB_AX_Y(z)*CROSSB_ZY;
            temp += d*g*laplnorm*interaction;
            kfd[1][i][edge_list[index]][j] += temp/3.;
            fdsub1[1][i][edge_list[index]][j] = fd[1][i][edge_list[index]][j] + temp;
          }
          if(i<=1 || i>=N-2) continue;
          for(index=0;index<4;index++)
          {
            // x=0,1,N-2,N-1
            temp = d*(-rescale_s*adsub2*fdsub2[1][edge_list[index]][i][j]/asub2 + laplnorm*laplb(1,edge_list[index],i,j,kloop_arg,FSUB));
            interaction  = DERIVB_AX_X(x)*CROSSB_XX + DERIVB_AX_X(y)*CROSSB_YX + DERIVB_AX_X(z)*CROSSB_ZX;
            temp += d*g*laplnorm*interaction;
            kfd[1][edge_list[index]][i][j] += temp/3.;
            fdsub1[1][edge_list[index]][i][j] = fd[1][edge_list[index]][i][j] + temp;
          }
        }
      break;
    }
    case 4:
    {
      laplnorm = 1./pow(asub1,2.*rescale_s);
      #pragma omp parallel for num_threads(num_thread) private(j,k,temp,interaction)
      for(i=2;i<N-2;i++)
        for(j=2;j<N-2;j++)
          for(k=2;k<N-2;k++)
          {
            temp = d*(-rescale_s*adsub1*GAUGEDSUB1(0)/asub1 + laplnorm*lapl(1,i,j,k,kloop_arg,FSUB));
            interaction  = DERIV_AX(x)*CROSSX + DERIV_AX(y)*CROSSY + DERIV_AX(z)*CROSSZ;
            temp += d*g*laplnorm*interaction;
            KGAUGED(0) += temp/6.;
          }
      #pragma omp parallel for num_threads(num_thread) private(j,temp,interaction,index)
      for(i=0;i<N;i++) // Index in x direction
        for(j=0;j<N;j++) // Index in y direction
        {
          for(index=0;index<4;index++)
          {
            // z=0,1,N-2,N-1
            temp = d*(-rescale_s*adsub1*fdsub1[1][i][j][edge_list[index]]/asub1 + laplnorm*laplb(1,i,j,edge_list[index],kloop_arg,FSUB));
            interaction  = DERIVB_AX_Z(x)*CROSSB_XZ + DERIVB_AX_Z(y)*CROSSB_YZ + DERIVB_AX_Z(z)*CROSSB_ZZ;
            temp += d*g*laplnorm*interaction;
            kfd[1][i][j][edge_list[index]] += temp/6.;
          }
          if(j<=1 || j>=N-2) continue;
          for(index=0;index<4;index++)
          {
            // y=0,1,N-2,N-1
            temp = d*(-rescale_s*adsub1*fdsub1[1][i][edge_list[index]][j]/asub1 + laplnorm*laplb(1,i,edge_list[index],j,kloop_arg,FSUB));
            interaction  = DERIVB_AX_Y(x)*CROSSB_XY + DERIVB_AX_Y(y)*CROSSB_YY + DERIVB_AX_Y(z)*CROSSB_ZY;
            temp += d*g*laplnorm*interaction;
            kfd[1][i][edge_list[index]][j] += temp/6.;
          }
          if(i<=1 || i>=N-2) continue;
          for(index=0;index<4;index++)
          {
            // x=0,1,N-2,N-1
            temp = d*(-rescale_s*adsub1*fdsub1[1][edge_list[index]][i][j]/asub1 + laplnorm*laplb(1,edge_list[index],i,j,kloop_arg,FSUB));
            interaction  = DERIVB_AX_X(x)*CROSSB_XX + DERIVB_AX_X(y)*CROSSB_YX + DERIVB_AX_X(z)*CROSSB_ZX;
            temp += d*g*laplnorm*interaction;
            kfd[1][edge_list[index]][i][j] += temp/6.;
          }
        }
      break;
    }
    default:
    {
      printf("Error: kloop takes unknown value\n");
      break;
    }
  }
}


// Evolve Time Derivative of Gauge Field (i component)
void evolve_derivs_gfi(double d, int kloop_arg, SUB_ARG, double kfd[nflds][N][N][N])
{
  int i=0,j=0,k=0,fld,fldg;
  double laplnorm,asnorm,asnorm2,temp,interaction;
  int edge_list[] = {0,1,N-2,N-1};
  int index;

  for(fld=2;fld<5;fld++)
  {
    fldg = fld - 1;
    switch(kloop_arg)
    {
      case 1:
      {
        laplnorm = 1./pow(a,2.*rescale_s);
        asnorm = 1./pow(a,rescale_s);
        asnorm2 = pow(a,rescale_s);
        #pragma omp parallel for num_threads(num_thread) private(j,k,temp,interaction)
        for(i=2;i<N-2;i++)
          for(j=2;j<N-2;j++)
            for(k=2;k<N-2;k++)
            {
              temp = d*(-rescale_s*ad*GAUGED(fldg)/a + laplnorm*lapl(fld,i,j,k,kloop_arg,FSUB));
              switch(fldg)
              {
                case 1:
                {
                  interaction  = asnorm*AXIOND*CROSSX;
                  interaction -= laplnorm*( DERIV_AX(y)*(asnorm2*GAUGED(3)-DERIV_GF(z,0)) - DERIV_AX(z)*(asnorm2*GAUGED(2)-DERIV_GF(y,0)) );
                  break;
                }
                case 2:
                {
                  interaction  = asnorm*AXIOND*CROSSY;
                  interaction -= laplnorm*( DERIV_AX(z)*(asnorm2*GAUGED(1)-DERIV_GF(x,0)) - DERIV_AX(x)*(asnorm2*GAUGED(3)-DERIV_GF(z,0)) );
                  break;
                }
                case 3:
                {
                  interaction  = asnorm*AXIOND*CROSSZ;
                  interaction -= laplnorm*( DERIV_AX(x)*(asnorm2*GAUGED(2)-DERIV_GF(y,0)) - DERIV_AX(y)*(asnorm2*GAUGED(1)-DERIV_GF(x,0)) );
                  break;
                }
                default:
                {
                  printf("Error: fldg takes unknown value\n");
                  break;
                }
              }
              temp += d*g*interaction;
              KGAUGED(fldg) = temp/6.;
              GAUGEDSUB1(fldg) = GAUGED(fldg) + temp/2.;
            }
        #pragma omp parallel for num_threads(num_thread) private(j,temp,interaction,index)
        for(i=0;i<N;i++) // Index in x direction
          for(j=0;j<N;j++) // Index in y direction
          {
            for(index=0;index<4;index++)
            {
              // z=0,1,N-2,N-1
              temp = d*(-rescale_s*ad*fd[fld][i][j][edge_list[index]]/a + laplnorm*laplb(fld,i,j,edge_list[index],kloop_arg,FSUB));
              switch(fldg)
              {
                case 1:
                {
                  interaction  = asnorm*fd[0][i][j][edge_list[index]]*CROSSB_XZ;
                  interaction -= laplnorm*( DERIVB_AX_Z(y)*(asnorm2*fd[4][i][j][edge_list[index]]-DERIVB_GF_Z(z,0)) - DERIVB_AX_Z(z)*(asnorm2*fd[3][i][j][edge_list[index]]-DERIVB_GF_Z(y,0)) );
                  break;
                }
                case 2:
                {
                  interaction  = asnorm*fd[0][i][j][edge_list[index]]*CROSSB_YZ;
                  interaction -= laplnorm*( DERIVB_AX_Z(z)*(asnorm2*fd[2][i][j][edge_list[index]]-DERIVB_GF_Z(x,0)) - DERIVB_AX_Z(x)*(asnorm2*fd[4][i][j][edge_list[index]]-DERIVB_GF_Z(z,0)) );
                  break;
                }
                case 3:
                {
                  interaction  = asnorm*fd[0][i][j][edge_list[index]]*CROSSB_ZZ;
                  interaction -= laplnorm*( DERIVB_AX_Z(x)*(asnorm2*fd[3][i][j][edge_list[index]]-DERIVB_GF_Z(y,0)) - DERIVB_AX_Z(y)*(asnorm2*fd[2][i][j][edge_list[index]]-DERIVB_GF_Z(x,0)) );
                  break;
                }
                default:
                {
                  printf("Error: fldg takes unknown value\n");
                  break;
                }
              }
              temp += d*g*interaction;
              kfd[fld][i][j][edge_list[index]] = temp/6.;
              fdsub1[fld][i][j][edge_list[index]] = fd[fld][i][j][edge_list[index]] + temp/2.;
            }
            if(j<=1 || j>=N-2) continue;
            for(index=0;index<4;index++)
            {
              // y=0,1,N-2,N-1
              temp = d*(-rescale_s*ad*fd[fld][i][edge_list[index]][j]/a + laplnorm*laplb(fld,i,edge_list[index],j,kloop_arg,FSUB));
              switch(fldg)
              {
                case 1:
                {
                  interaction  = asnorm*fd[0][i][edge_list[index]][j]*CROSSB_XY;
                  interaction -= laplnorm*( DERIVB_AX_Y(y)*(asnorm2*fd[4][i][edge_list[index]][j]-DERIVB_GF_Y(z,0)) - DERIVB_AX_Y(z)*(asnorm2*fd[3][i][edge_list[index]][j]-DERIVB_GF_Y(y,0)) );
                  break;
                }
                case 2:
                {
                  interaction  = asnorm*fd[0][i][edge_list[index]][j]*CROSSB_YY;
                  interaction -= laplnorm*( DERIVB_AX_Y(z)*(asnorm2*fd[2][i][edge_list[index]][j]-DERIVB_GF_Y(x,0)) - DERIVB_AX_Y(x)*(asnorm2*fd[4][i][edge_list[index]][j]-DERIVB_GF_Y(z,0)) );
                  break;
                }
                case 3:
                {
                  interaction  = asnorm*fd[0][i][edge_list[index]][j]*CROSSB_ZY;
                  interaction -= laplnorm*( DERIVB_AX_Y(x)*(asnorm2*fd[3][i][edge_list[index]][j]-DERIVB_GF_Y(y,0)) - DERIVB_AX_Y(y)*(asnorm2*fd[2][i][edge_list[index]][j]-DERIVB_GF_Y(x,0)) );
                  break;
                }
                default:
                {
                  printf("Error: fldg takes unknown value\n");
                  break;
                }
              }
              temp += d*g*interaction;
              kfd[fld][i][edge_list[index]][j] = temp/6.;
              fdsub1[fld][i][edge_list[index]][j] = fd[fld][i][edge_list[index]][j] + temp/2.;
            }
            if(i<=1 || i>=N-2) continue;
            for(index=0;index<4;index++)
            {
              // x=0,1,N-2,N-1
              temp = d*(-rescale_s*ad*fd[fld][edge_list[index]][i][j]/a + laplnorm*laplb(fld,edge_list[index],i,j,kloop_arg,FSUB));
              switch(fldg)
              {
                case 1:
                {
                  interaction  = asnorm*fd[0][edge_list[index]][i][j]*CROSSB_XX;
                  interaction -= laplnorm*( DERIVB_AX_X(y)*(asnorm2*fd[4][edge_list[index]][i][j]-DERIVB_GF_X(z,0)) - DERIVB_AX_X(z)*(asnorm2*fd[3][edge_list[index]][i][j]-DERIVB_GF_X(y,0)) );
                  break;
                }
                case 2:
                {
                  interaction  = asnorm*fd[0][edge_list[index]][i][j]*CROSSB_YX;
                  interaction -= laplnorm*( DERIVB_AX_X(z)*(asnorm2*fd[2][edge_list[index]][i][j]-DERIVB_GF_X(x,0)) - DERIVB_AX_X(x)*(asnorm2*fd[4][edge_list[index]][i][j]-DERIVB_GF_X(z,0)) );
                  break;
                }
                case 3:
                {
                  interaction  = asnorm*fd[0][edge_list[index]][i][j]*CROSSB_ZX;
                  interaction -= laplnorm*( DERIVB_AX_X(x)*(asnorm2*fd[3][edge_list[index]][i][j]-DERIVB_GF_X(y,0)) - DERIVB_AX_X(y)*(asnorm2*fd[2][edge_list[index]][i][j]-DERIVB_GF_X(x,0)) );
                  break;
                }
                default:
                {
                  printf("Error: fldg takes unknown value\n");
                  break;
                }
              }
              temp += d*g*interaction;
              kfd[fld][edge_list[index]][i][j] = temp/6.;
              fdsub1[fld][edge_list[index]][i][j] = fd[fld][edge_list[index]][i][j] + temp/2.;
            }
          }
        break;
      }
      case 2:
      {
        laplnorm = 1./pow(asub1,2.*rescale_s);
        asnorm = 1./pow(asub1,rescale_s);
        asnorm2 = pow(asub1,rescale_s);
        #pragma omp parallel for num_threads(num_thread) private(j,k,temp,interaction)
        for(i=2;i<N-2;i++)
          for(j=2;j<N-2;j++)
            for(k=2;k<N-2;k++)
            {
              temp = d*(-rescale_s*adsub1*GAUGEDSUB1(fldg)/asub1 + laplnorm*lapl(fld,i,j,k,kloop_arg,FSUB));
              switch(fldg)
              {
                case 1:
                {
                  interaction  = asnorm*AXIONDSUB1*CROSSX;
                  interaction -= laplnorm*( DERIV_AX(y)*(asnorm2*GAUGEDSUB1(3)-DERIV_GF(z,0)) - DERIV_AX(z)*(asnorm2*GAUGEDSUB1(2)-DERIV_GF(y,0)) );
                  break;
                }
                case 2:
                {
                  interaction  = asnorm*AXIONDSUB1*CROSSY;
                  interaction -= laplnorm*( DERIV_AX(z)*(asnorm2*GAUGEDSUB1(1)-DERIV_GF(x,0)) - DERIV_AX(x)*(asnorm2*GAUGEDSUB1(3)-DERIV_GF(z,0)) );
                  break;
                }
                case 3:
                {
                  interaction  = asnorm*AXIONDSUB1*CROSSZ;
                  interaction -= laplnorm*( DERIV_AX(x)*(asnorm2*GAUGEDSUB1(2)-DERIV_GF(y,0)) - DERIV_AX(y)*(asnorm2*GAUGEDSUB1(1)-DERIV_GF(x,0)) );
                  break;
                }
                default:
                {
                  printf("Error: fldg takes unknown value\n");
                  break;
                }
              }
              temp += d*g*interaction;
              KGAUGED(fldg) += temp/3.;
              GAUGEDSUB2(fldg) = GAUGED(fldg) + temp/2.;
            }
        #pragma omp parallel for num_threads(num_thread) private(j,temp,interaction,index)
        for(i=0;i<N;i++) // Index in x direction
          for(j=0;j<N;j++) // Index in y direction
          {
            for(index=0;index<4;index++)
            {
              // z=0,1,N-2,N-1
              temp = d*(-rescale_s*adsub1*fdsub1[fld][i][j][edge_list[index]]/asub1 + laplnorm*laplb(fld,i,j,edge_list[index],kloop_arg,FSUB));
              switch(fldg)
              {
                case 1:
                {
                  interaction  = asnorm*fdsub1[0][i][j][edge_list[index]]*CROSSB_XZ;
                  interaction -= laplnorm*( DERIVB_AX_Z(y)*(asnorm2*fdsub1[4][i][j][edge_list[index]]-DERIVB_GF_Z(z,0)) - DERIVB_AX_Z(z)*(asnorm2*fdsub1[3][i][j][edge_list[index]]-DERIVB_GF_Z(y,0)) );
                  break;
                }
                case 2:
                {
                  interaction  = asnorm*fdsub1[0][i][j][edge_list[index]]*CROSSB_YZ;
                  interaction -= laplnorm*( DERIVB_AX_Z(z)*(asnorm2*fdsub1[2][i][j][edge_list[index]]-DERIVB_GF_Z(x,0)) - DERIVB_AX_Z(x)*(asnorm2*fdsub1[4][i][j][edge_list[index]]-DERIVB_GF_Z(z,0)) );
                  break;
                }
                case 3:
                {
                  interaction  = asnorm*fdsub1[0][i][j][edge_list[index]]*CROSSB_ZZ;
                  interaction -= laplnorm*( DERIVB_AX_Z(x)*(asnorm2*fdsub1[3][i][j][edge_list[index]]-DERIVB_GF_Z(y,0)) - DERIVB_AX_Z(y)*(asnorm2*fdsub1[2][i][j][edge_list[index]]-DERIVB_GF_Z(x,0)) );
                  break;
                }
                default:
                {
                  printf("Error: fldg takes unknown value\n");
                  break;
                }
              }
              temp += d*g*interaction;
              kfd[fld][i][j][edge_list[index]] += temp/3.;
              fdsub2[fld][i][j][edge_list[index]] = fd[fld][i][j][edge_list[index]] + temp/2.;
            }
            if(j<=1 || j>=N-2) continue;
            for(index=0;index<4;index++)
            {
              // y=0,1,N-2,N-1
              temp = d*(-rescale_s*adsub1*fdsub1[fld][i][edge_list[index]][j]/asub1 + laplnorm*laplb(fld,i,edge_list[index],j,kloop_arg,FSUB));
              switch(fldg)
              {
                case 1:
                {
                  interaction  = asnorm*fdsub1[0][i][edge_list[index]][j]*CROSSB_XY;
                  interaction -= laplnorm*( DERIVB_AX_Y(y)*(asnorm2*fdsub1[4][i][edge_list[index]][j]-DERIVB_GF_Y(z,0)) - DERIVB_AX_Y(z)*(asnorm2*fdsub1[3][i][edge_list[index]][j]-DERIVB_GF_Y(y,0)) );
                  break;
                }
                case 2:
                {
                  interaction  = asnorm*fdsub1[0][i][edge_list[index]][j]*CROSSB_YY;
                  interaction -= laplnorm*( DERIVB_AX_Y(z)*(asnorm2*fdsub1[2][i][edge_list[index]][j]-DERIVB_GF_Y(x,0)) - DERIVB_AX_Y(x)*(asnorm2*fdsub1[4][i][edge_list[index]][j]-DERIVB_GF_Y(z,0)) );
                  break;
                }
                case 3:
                {
                  interaction  = asnorm*fdsub1[0][i][edge_list[index]][j]*CROSSB_ZY;
                  interaction -= laplnorm*( DERIVB_AX_Y(x)*(asnorm2*fdsub1[3][i][edge_list[index]][j]-DERIVB_GF_Y(y,0)) - DERIVB_AX_Y(y)*(asnorm2*fdsub1[2][i][edge_list[index]][j]-DERIVB_GF_Y(x,0)) );
                  break;
                }
                default:
                {
                  printf("Error: fldg takes unknown value\n");
                  break;
                }
              }
              temp += d*g*interaction;
              kfd[fld][i][edge_list[index]][j] += temp/3.;
              fdsub2[fld][i][edge_list[index]][j] = fd[fld][i][edge_list[index]][j] + temp/2.;
            }
            if(i<=1 || i>=N-2) continue;
            for(index=0;index<4;index++)
            {
              // x=0,1,N-2,N-1
              temp = d*(-rescale_s*adsub1*fdsub1[fld][edge_list[index]][i][j]/asub1 + laplnorm*laplb(fld,edge_list[index],i,j,kloop_arg,FSUB));
              switch(fldg)
              {
                case 1:
                {
                  interaction  = asnorm*fdsub1[0][edge_list[index]][i][j]*CROSSB_XX;
                  interaction -= laplnorm*( DERIVB_AX_X(y)*(asnorm2*fdsub1[4][edge_list[index]][i][j]-DERIVB_GF_X(z,0)) - DERIVB_AX_X(z)*(asnorm2*fdsub1[3][edge_list[index]][i][j]-DERIVB_GF_X(y,0)) );
                  break;
                }
                case 2:
                {
                  interaction  = asnorm*fdsub1[0][edge_list[index]][i][j]*CROSSB_YX;
                  interaction -= laplnorm*( DERIVB_AX_X(z)*(asnorm2*fdsub1[2][edge_list[index]][i][j]-DERIVB_GF_X(x,0)) - DERIVB_AX_X(x)*(asnorm2*fdsub1[4][edge_list[index]][i][j]-DERIVB_GF_X(z,0)) );
                  break;
                }
                case 3:
                {
                  interaction  = asnorm*fdsub1[0][edge_list[index]][i][j]*CROSSB_ZX;
                  interaction -= laplnorm*( DERIVB_AX_X(x)*(asnorm2*fdsub1[3][edge_list[index]][i][j]-DERIVB_GF_X(y,0)) - DERIVB_AX_X(y)*(asnorm2*fdsub1[2][edge_list[index]][i][j]-DERIVB_GF_X(x,0)) );
                  break;
                }
                default:
                {
                  printf("Error: fldg takes unknown value\n");
                  break;
                }
              }
              temp += d*g*interaction;
              kfd[fld][edge_list[index]][i][j] += temp/3.;
              fdsub2[fld][edge_list[index]][i][j] = fd[fld][edge_list[index]][i][j] + temp/2.;
            }
          }
        break;
      }
      case 3:
      {
        laplnorm = 1./pow(asub2,2.*rescale_s);
        asnorm = 1./pow(asub2,rescale_s);
        asnorm2 = pow(asub2,rescale_s);
        #pragma omp parallel for num_threads(num_thread) private(j,k,temp,interaction)
        for(i=2;i<N-2;i++)
          for(j=2;j<N-2;j++)
            for(k=2;k<N-2;k++)
            {
              temp = d*(-rescale_s*adsub2*GAUGEDSUB2(fldg)/asub2 + laplnorm*lapl(fld,i,j,k,kloop_arg,FSUB));
              switch(fldg)
              {
                case 1:
                {
                  interaction  = asnorm*AXIONDSUB2*CROSSX;
                  interaction -= laplnorm*( DERIV_AX(y)*(asnorm2*GAUGEDSUB2(3)-DERIV_GF(z,0)) - DERIV_AX(z)*(asnorm2*GAUGEDSUB2(2)-DERIV_GF(y,0)) );
                  break;
                }
                case 2:
                {
                  interaction  = asnorm*AXIONDSUB2*CROSSY;
                  interaction -= laplnorm*( DERIV_AX(z)*(asnorm2*GAUGEDSUB2(1)-DERIV_GF(x,0)) - DERIV_AX(x)*(asnorm2*GAUGEDSUB2(3)-DERIV_GF(z,0)) );
                  break;
                }
                case 3:
                {
                  interaction  = asnorm*AXIONDSUB2*CROSSZ;
                  interaction -= laplnorm*( DERIV_AX(x)*(asnorm2*GAUGEDSUB2(2)-DERIV_GF(y,0)) - DERIV_AX(y)*(asnorm2*GAUGEDSUB2(1)-DERIV_GF(x,0)) );
                  break;
                }
                default:
                {
                  printf("Error: fldg takes unknown value\n");
                  break;
                }
              }
              temp += d*g*interaction;
              KGAUGED(fldg) += temp/3.;
              GAUGEDSUB1(fldg) = GAUGED(fldg) + temp;
            }
        #pragma omp parallel for num_threads(num_thread) private(j,temp,interaction,index)
        for(i=0;i<N;i++) // Index in x direction
          for(j=0;j<N;j++) // Index in y direction
          {
            for(index=0;index<4;index++)
            {
              // z=0,1,N-2,N-1
              temp = d*(-rescale_s*adsub2*fdsub2[fld][i][j][edge_list[index]]/asub2 + laplnorm*laplb(fld,i,j,edge_list[index],kloop_arg,FSUB));
              switch(fldg)
              {
                case 1:
                {
                  interaction  = asnorm*fdsub2[0][i][j][edge_list[index]]*CROSSB_XZ;
                  interaction -= laplnorm*( DERIVB_AX_Z(y)*(asnorm2*fdsub2[4][i][j][edge_list[index]]-DERIVB_GF_Z(z,0)) - DERIVB_AX_Z(z)*(asnorm2*fdsub2[3][i][j][edge_list[index]]-DERIVB_GF_Z(y,0)) );
                  break;
                }
                case 2:
                {
                  interaction  = asnorm*fdsub2[0][i][j][edge_list[index]]*CROSSB_YZ;
                  interaction -= laplnorm*( DERIVB_AX_Z(z)*(asnorm2*fdsub2[2][i][j][edge_list[index]]-DERIVB_GF_Z(x,0)) - DERIVB_AX_Z(x)*(asnorm2*fdsub2[4][i][j][edge_list[index]]-DERIVB_GF_Z(z,0)) );
                  break;
                }
                case 3:
                {
                  interaction  = asnorm*fdsub2[0][i][j][edge_list[index]]*CROSSB_ZZ;
                  interaction -= laplnorm*( DERIVB_AX_Z(x)*(asnorm2*fdsub2[3][i][j][edge_list[index]]-DERIVB_GF_Z(y,0)) - DERIVB_AX_Z(y)*(asnorm2*fdsub2[2][i][j][edge_list[index]]-DERIVB_GF_Z(x,0)) );
                  break;
                }
                default:
                {
                  printf("Error: fldg takes unknown value\n");
                  break;
                }
              }
              temp += d*g*interaction;
              kfd[fld][i][j][edge_list[index]] += temp/3.;
              fdsub1[fld][i][j][edge_list[index]] = fd[fld][i][j][edge_list[index]] + temp;
            }
            if(j<=1 || j>=N-2) continue;
            for(index=0;index<4;index++)
            {
              // y=0,1,N-2,N-1
              temp = d*(-rescale_s*adsub2*fdsub2[fld][i][edge_list[index]][j]/asub2 + laplnorm*laplb(fld,i,edge_list[index],j,kloop_arg,FSUB));
              switch(fldg)
              {
                case 1:
                {
                  interaction  = asnorm*fdsub2[0][i][edge_list[index]][j]*CROSSB_XY;
                  interaction -= laplnorm*( DERIVB_AX_Y(y)*(asnorm2*fdsub2[4][i][edge_list[index]][j]-DERIVB_GF_Y(z,0)) - DERIVB_AX_Y(z)*(asnorm2*fdsub2[3][i][edge_list[index]][j]-DERIVB_GF_Y(y,0)) );
                  break;
                }
                case 2:
                {
                  interaction  = asnorm*fdsub2[0][i][edge_list[index]][j]*CROSSB_YY;
                  interaction -= laplnorm*( DERIVB_AX_Y(z)*(asnorm2*fdsub2[2][i][edge_list[index]][j]-DERIVB_GF_Y(x,0)) - DERIVB_AX_Y(x)*(asnorm2*fdsub2[4][i][edge_list[index]][j]-DERIVB_GF_Y(z,0)) );
                  break;
                }
                case 3:
                {
                  interaction  = asnorm*fdsub2[0][i][edge_list[index]][j]*CROSSB_ZY;
                  interaction -= laplnorm*( DERIVB_AX_Y(x)*(asnorm2*fdsub2[3][i][edge_list[index]][j]-DERIVB_GF_Y(y,0)) - DERIVB_AX_Y(y)*(asnorm2*fdsub2[2][i][edge_list[index]][j]-DERIVB_GF_Y(x,0)) );
                  break;
                }
                default:
                {
                  printf("Error: fldg takes unknown value\n");
                  break;
                }
              }
              temp += d*g*interaction;
              kfd[fld][i][edge_list[index]][j] += temp/3.;
              fdsub1[fld][i][edge_list[index]][j] = fd[fld][i][edge_list[index]][j] + temp;
            }
            if(i<=1 || i>=N-2) continue;
            for(index=0;index<4;index++)
            {
              // x=0,1,N-2,N-1
              temp = d*(-rescale_s*adsub2*fdsub2[fld][edge_list[index]][i][j]/asub2 + laplnorm*laplb(fld,edge_list[index],i,j,kloop_arg,FSUB));
              switch(fldg)
              {
                case 1:
                {
                  interaction  = asnorm*fdsub2[0][edge_list[index]][i][j]*CROSSB_XX;
                  interaction -= laplnorm*( DERIVB_AX_X(y)*(asnorm2*fdsub2[4][edge_list[index]][i][j]-DERIVB_GF_X(z,0)) - DERIVB_AX_X(z)*(asnorm2*fdsub2[3][edge_list[index]][i][j]-DERIVB_GF_X(y,0)) );
                  break;
                }
                case 2:
                {
                  interaction  = asnorm*fdsub2[0][edge_list[index]][i][j]*CROSSB_YX;
                  interaction -= laplnorm*( DERIVB_AX_X(z)*(asnorm2*fdsub2[2][edge_list[index]][i][j]-DERIVB_GF_X(x,0)) - DERIVB_AX_X(x)*(asnorm2*fdsub2[4][edge_list[index]][i][j]-DERIVB_GF_X(z,0)) );
                  break;
                }
                case 3:
                {
                  interaction  = asnorm*fdsub2[0][edge_list[index]][i][j]*CROSSB_ZX;
                  interaction -= laplnorm*( DERIVB_AX_X(x)*(asnorm2*fdsub2[3][edge_list[index]][i][j]-DERIVB_GF_X(y,0)) - DERIVB_AX_X(y)*(asnorm2*fdsub2[2][edge_list[index]][i][j]-DERIVB_GF_X(x,0)) );
                  break;
                }
                default:
                {
                  printf("Error: fldg takes unknown value\n");
                  break;
                }
              }
              temp += d*g*interaction;
              kfd[fld][edge_list[index]][i][j] += temp/3.;
              fdsub1[fld][edge_list[index]][i][j] = fd[fld][edge_list[index]][i][j] + temp;
            }
          }
        break;
      }
      case 4:
      {
        laplnorm = 1./pow(asub1,2.*rescale_s);
        asnorm = 1./pow(asub1,rescale_s);
        asnorm2 = pow(asub1,rescale_s);
        #pragma omp parallel for num_threads(num_thread) private(j,k,temp,interaction)
        for(i=2;i<N-2;i++)
          for(j=2;j<N-2;j++)
            for(k=2;k<N-2;k++)
            {
              temp = d*(-rescale_s*adsub1*GAUGEDSUB1(fldg)/asub1 + laplnorm*lapl(fld,i,j,k,kloop_arg,FSUB));
              switch(fldg)
              {
                case 1:
                {
                  interaction  = asnorm*AXIONDSUB1*CROSSX;
                  interaction -= laplnorm*( DERIV_AX(y)*(asnorm2*GAUGEDSUB1(3)-DERIV_GF(z,0)) - DERIV_AX(z)*(asnorm2*GAUGEDSUB1(2)-DERIV_GF(y,0)) );
                  break;
                }
                case 2:
                {
                  interaction  = asnorm*AXIONDSUB1*CROSSY;
                  interaction -= laplnorm*( DERIV_AX(z)*(asnorm2*GAUGEDSUB1(1)-DERIV_GF(x,0)) - DERIV_AX(x)*(asnorm2*GAUGEDSUB1(3)-DERIV_GF(z,0)) );
                  break;
                }
                case 3:
                {
                  interaction  = asnorm*AXIONDSUB1*CROSSZ;
                  interaction -= laplnorm*( DERIV_AX(x)*(asnorm2*GAUGEDSUB1(2)-DERIV_GF(y,0)) - DERIV_AX(y)*(asnorm2*GAUGEDSUB1(1)-DERIV_GF(x,0)) );
                  break;
                }
                default:
                {
                  printf("Error: fldg takes unknown value\n");
                  break;
                }
              }
              temp += d*g*interaction;
              KGAUGED(fldg) += temp/6.;
            }
        #pragma omp parallel for num_threads(num_thread) private(j,temp,interaction,index)
        for(i=0;i<N;i++) // Index in x direction
          for(j=0;j<N;j++) // Index in y direction
          {
            for(index=0;index<4;index++)
            {
              // z=0,1,N-2,N-1
              temp = d*(-rescale_s*adsub1*fdsub1[fld][i][j][edge_list[index]]/asub1 + laplnorm*laplb(fld,i,j,edge_list[index],kloop_arg,FSUB));
              switch(fldg)
              {
                case 1:
                {
                  interaction  = asnorm*fdsub1[0][i][j][edge_list[index]]*CROSSB_XZ;
                  interaction -= laplnorm*( DERIVB_AX_Z(y)*(asnorm2*fdsub1[4][i][j][edge_list[index]]-DERIVB_GF_Z(z,0)) - DERIVB_AX_Z(z)*(asnorm2*fdsub1[3][i][j][edge_list[index]]-DERIVB_GF_Z(y,0)) );
                  break;
                }
                case 2:
                {
                  interaction  = asnorm*fdsub1[0][i][j][edge_list[index]]*CROSSB_YZ;
                  interaction -= laplnorm*( DERIVB_AX_Z(z)*(asnorm2*fdsub1[2][i][j][edge_list[index]]-DERIVB_GF_Z(x,0)) - DERIVB_AX_Z(x)*(asnorm2*fdsub1[4][i][j][edge_list[index]]-DERIVB_GF_Z(z,0)) );
                  break;
                }
                case 3:
                {
                  interaction  = asnorm*fdsub1[0][i][j][edge_list[index]]*CROSSB_ZZ;
                  interaction -= laplnorm*( DERIVB_AX_Z(x)*(asnorm2*fdsub1[3][i][j][edge_list[index]]-DERIVB_GF_Z(y,0)) - DERIVB_AX_Z(y)*(asnorm2*fdsub1[2][i][j][edge_list[index]]-DERIVB_GF_Z(x,0)) );
                  break;
                }
                default:
                {
                  printf("Error: fldg takes unknown value\n");
                  break;
                }
              }
              temp += d*g*interaction;
              kfd[fld][i][j][edge_list[index]] += temp/6.;
            }
            if(j<=1 || j>=N-2) continue;
            for(index=0;index<4;index++)
            {
              // y=0,1,N-2,N-1
              temp = d*(-rescale_s*adsub1*fdsub1[fld][i][edge_list[index]][j]/asub1 + laplnorm*laplb(fld,i,edge_list[index],j,kloop_arg,FSUB));
              switch(fldg)
              {
                case 1:
                {
                  interaction  = asnorm*fdsub1[0][i][edge_list[index]][j]*CROSSB_XY;
                  interaction -= laplnorm*( DERIVB_AX_Y(y)*(asnorm2*fdsub1[4][i][edge_list[index]][j]-DERIVB_GF_Y(z,0)) - DERIVB_AX_Y(z)*(asnorm2*fdsub1[3][i][edge_list[index]][j]-DERIVB_GF_Y(y,0)) );
                  break;
                }
                case 2:
                {
                  interaction  = asnorm*fdsub1[0][i][edge_list[index]][j]*CROSSB_YY;
                  interaction -= laplnorm*( DERIVB_AX_Y(z)*(asnorm2*fdsub1[2][i][edge_list[index]][j]-DERIVB_GF_Y(x,0)) - DERIVB_AX_Y(x)*(asnorm2*fdsub1[4][i][edge_list[index]][j]-DERIVB_GF_Y(z,0)) );
                  break;
                }
                case 3:
                {
                  interaction  = asnorm*fdsub1[0][i][edge_list[index]][j]*CROSSB_ZY;
                  interaction -= laplnorm*( DERIVB_AX_Y(x)*(asnorm2*fdsub1[3][i][edge_list[index]][j]-DERIVB_GF_Y(y,0)) - DERIVB_AX_Y(y)*(asnorm2*fdsub1[2][i][edge_list[index]][j]-DERIVB_GF_Y(x,0)) );
                  break;
                }
                default:
                {
                  printf("Error: fldg takes unknown value\n");
                  break;
                }
              }
              temp += d*g*interaction;
              kfd[fld][i][edge_list[index]][j] += temp/6.;
            }
            if(i<=1 || i>=N-2) continue;
            for(index=0;index<4;index++)
            {
              // x=0,1,N-2,N-1
              temp = d*(-rescale_s*adsub1*fdsub1[fld][edge_list[index]][i][j]/asub1 + laplnorm*laplb(fld,edge_list[index],i,j,kloop_arg,FSUB));
              switch(fldg)
              {
                case 1:
                {
                  interaction  = asnorm*fdsub1[0][edge_list[index]][i][j]*CROSSB_XX;
                  interaction -= laplnorm*( DERIVB_AX_X(y)*(asnorm2*fdsub1[4][edge_list[index]][i][j]-DERIVB_GF_X(z,0)) - DERIVB_AX_X(z)*(asnorm2*fdsub1[3][edge_list[index]][i][j]-DERIVB_GF_X(y,0)) );
                  break;
                }
                case 2:
                {
                  interaction  = asnorm*fdsub1[0][edge_list[index]][i][j]*CROSSB_YX;
                  interaction -= laplnorm*( DERIVB_AX_X(z)*(asnorm2*fdsub1[2][edge_list[index]][i][j]-DERIVB_GF_X(x,0)) - DERIVB_AX_X(x)*(asnorm2*fdsub1[4][edge_list[index]][i][j]-DERIVB_GF_X(z,0)) );
                  break;
                }
                case 3:
                {
                  interaction  = asnorm*fdsub1[0][edge_list[index]][i][j]*CROSSB_ZX;
                  interaction -= laplnorm*( DERIVB_AX_X(x)*(asnorm2*fdsub1[3][edge_list[index]][i][j]-DERIVB_GF_X(y,0)) - DERIVB_AX_X(y)*(asnorm2*fdsub1[2][edge_list[index]][i][j]-DERIVB_GF_X(x,0)) );
                  break;
                }
                default:
                {
                  printf("Error: fldg takes unknown value\n");
                  break;
                }
              }
              temp += d*g*interaction;
              kfd[fld][edge_list[index]][i][j] += temp/6.;
            }
          }
        break;
      }
      default:
      {
        printf("Error: kloop takes unknown value\n");
        break;
      }
    }
  }
}


// Evolve Scale Factor
void evolve_scale(double d, int kloop_arg)
{
  switch(kloop_arg)
  {
    case 1:
    {
      ka = d*ad/6.;
      asub1 = a + d*ad/2.;
      break;
    }
    case 2:
    {
      ka += d*adsub1/3.;
      asub2 = a + d*adsub1/2.;
      break;
    }
    case 3:
    {
      ka += d*adsub2/3.;
      asub1 = a + d*adsub2;
      break;
    }
    case 4:
    {
      ka += d*adsub1/6.;
      break;
    }
    default:
    {
      printf("Error: kloop takes unknown value\n");
      break;
    }
  }
}


// Evolve Time Derivative of Scale Factor
void evolve_scaled(double d, int kloop_arg, SUB_ARG)
{
  double tempad;

  switch(kloop_arg)
  {
    case 1:
    {
      tempad = d*(-rescale_s*pw2(ad)/a + pow(a,3.)*(stress_3pressure(kloop_arg,SUB))/6.);
      kad = tempad/6.;
      adsub1 = ad + tempad/2.;
      break;
    }
    case 2:
    {
      tempad = d*(-rescale_s*pw2(adsub1)/asub1 + pow(asub1,3.)*(stress_3pressure(kloop_arg,SUB))/6.);
      kad += tempad/3.;
      adsub2 = ad + tempad/2.;
      break;
    }
    case 3:
    {
      tempad = d*(-rescale_s*pw2(adsub2)/asub2 + pow(asub2,3.)*(stress_3pressure(kloop_arg,SUB))/6.);
      kad += tempad/3.;
      adsub1 = ad + tempad;
      break;
    }
    case 4:
    {
      kad += d*(-rescale_s*pw2(adsub1)/asub1 + pow(asub1,3.)*(stress_3pressure(kloop_arg,SUB))/6.)/6.;
      break;
    }
    default:
    {
      printf("Error: kloop takes unknown value\n");
      break;
    }
  }
}



// Important Function of Evolve 1 Step
void evolve_rk4(double dt, SUB_ARG, KSUB_ARG)
{
  DECLARE_INDICES
  int kloop,fld;

  /* calculate delta */
  for(kloop=1;kloop<=4;kloop++)
  {
    // sub2->sub1を作るときにsub1を初期化しておく
    if(kloop==3)
    {
      for(int n=0;n<nflds;n++)
      for(int i=0;i<N;i++)
      for(int j=0;j<N;j++)
      for(int k=0;k<N;k++)
      {
        fsub1[n][i][j][k] = 0.;
        fdsub1[n][i][j][k] = 0.;
        asub1 = 0.;
        adsub1 = 0.;
      }
    }

    evolve_fields(dt,kloop,SUB,kf);
    evolve_derivs_ax(dt,kloop,SUB,kfd);
    evolve_derivs_gf0(dt,kloop,SUB,kfd);
    evolve_derivs_gfi(dt,kloop,SUB,kfd);
    evolve_scale (dt,kloop);
    evolve_scaled(dt,kloop,SUB);
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
