#ifndef _SPECTRA_
#define _SPECTRA_

#include "latticeeasy.hpp"

char name_[500]; //Filenames


/* Function calculating gauge amplitude */
/*----------------------------------------------------------------------*/
// |A+|^2 for k<N/2
double amplitude_gfposi(double px, double py, double pz, int i, int j, int k)
{
  double real,imag;
  double amp;

  if(px==0. && py==0.)
  {
    amp = (pw2(f[1][i][j][2*k]+f[2][i][j][2*k+1]) + pw2(f[1][i][j][2*k+1]-f[2][i][j][2*k]))/2.;
  }
  else
  {
    real = 0.;
    imag = 0.;

    // x-component
    real  = pola_x_re(px,py,pz)*f[1][i][j][2*k]   + pola_x_im(px,py,pz)*f[1][i][j][2*k+1];
    imag  = pola_x_re(px,py,pz)*f[1][i][j][2*k+1] - pola_x_im(px,py,pz)*f[1][i][j][2*k];

    // y-component
    real  += pola_y_re(px,py,pz)*f[2][i][j][2*k]   + pola_y_im(px,py,pz)*f[2][i][j][2*k+1];
    imag  += pola_y_re(px,py,pz)*f[2][i][j][2*k+1] - pola_y_im(px,py,pz)*f[2][i][j][2*k];

    // z-component
    real  += pola_z_re(px,py,pz)*f[3][i][j][2*k]   + pola_z_im(px,py,pz)*f[3][i][j][2*k+1];
    imag  += pola_z_re(px,py,pz)*f[3][i][j][2*k+1] - pola_z_im(px,py,pz)*f[3][i][j][2*k];

    // Add 3 directions
    amp = pw2(real) + pw2(imag);
  }

  return amp;
}


// |A+|^2 for k=N/2
double amplitude_gfposi_nyquist(double px, double py, int i, int j)
{
  double pz = (double)N/2.;

  double real,imag;
  double amp;

  if(px==0. && py==0.)
  {
    amp = (pw2(fnyquist_gf[0][i][2*j]+fnyquist_gf[1][i][2*j+1]) + pw2(fnyquist_gf[0][i][2*j+1]-fnyquist_gf[1][i][2*j]))/2.;
  }
  else
  {
    real = 0.;
    imag = 0.;

    // x-component
    real = pola_x_re(px,py,pz)*fnyquist_gf[0][i][2*j]   + pola_x_im(px,py,pz)*fnyquist_gf[0][i][2*j+1];
    imag = pola_x_re(px,py,pz)*fnyquist_gf[0][i][2*j+1] - pola_x_im(px,py,pz)*fnyquist_gf[0][i][2*j];

    // y-component
    real += pola_y_re(px,py,pz)*fnyquist_gf[1][i][2*j]   + pola_y_im(px,py,pz)*fnyquist_gf[1][i][2*j+1];
    imag += pola_y_re(px,py,pz)*fnyquist_gf[1][i][2*j+1] - pola_y_im(px,py,pz)*fnyquist_gf[1][i][2*j];

    // z-component
    real += pola_z_re(px,py,pz)*fnyquist_gf[2][i][2*j]   + pola_z_im(px,py,pz)*fnyquist_gf[2][i][2*j+1];
    imag += pola_z_re(px,py,pz)*fnyquist_gf[2][i][2*j+1] - pola_z_im(px,py,pz)*fnyquist_gf[2][i][2*j];

    // Add 3 directions
    amp = pw2(real) + pw2(imag);
  }

  return amp;
}


// |A-|^2 for k<N/2
double amplitude_gfnega(double px, double py, double pz, int i, int j, int k)
{
  double real,imag;
  double amp;

  if(px==0. && py==0.)
  {
    amp = (pw2(f[1][i][j][2*k]-f[2][i][j][2*k+1]) + pw2(f[1][i][j][2*k+1]+f[2][i][j][2*k]))/2.;
  }
  else
  {
    real = 0.;
    imag = 0.;

    // x-component
    real  = pola_x_re(px,py,pz)*f[1][i][j][2*k]   - pola_x_im(px,py,pz)*f[1][i][j][2*k+1];
    imag  = pola_x_re(px,py,pz)*f[1][i][j][2*k+1] + pola_x_im(px,py,pz)*f[1][i][j][2*k];

    // y-component
    real  += pola_y_re(px,py,pz)*f[2][i][j][2*k]   - pola_y_im(px,py,pz)*f[2][i][j][2*k+1];
    imag  += pola_y_re(px,py,pz)*f[2][i][j][2*k+1] + pola_y_im(px,py,pz)*f[2][i][j][2*k];

    // z-component
    real  += pola_z_re(px,py,pz)*f[3][i][j][2*k]   - pola_z_im(px,py,pz)*f[3][i][j][2*k+1];
    imag  += pola_z_re(px,py,pz)*f[3][i][j][2*k+1] + pola_z_im(px,py,pz)*f[3][i][j][2*k];

    // Add 3 directions
    amp = pw2(real) + pw2(imag);
  }

  return amp;
}


// |A-|^2 for k=N/2
double amplitude_gfnega_nyquist(double px, double py, int i, int j)
{
  double pz = (double)N/2.;

  double real,imag;
  double amp;

  if(px==0. && py==0.)
  {
    amp = (pw2(fnyquist_gf[0][i][2*j]-fnyquist_gf[1][i][2*j+1]) + pw2(fnyquist_gf[0][i][2*j+1]+fnyquist_gf[1][i][2*j]))/2.;
  }
  else
  {
    real = 0.;
    imag = 0.;

    // x-component
    real = pola_x_re(px,py,pz)*fnyquist_gf[0][i][2*j]   - pola_x_im(px,py,pz)*fnyquist_gf[0][i][2*j+1];
    imag = pola_x_re(px,py,pz)*fnyquist_gf[0][i][2*j+1] + pola_x_im(px,py,pz)*fnyquist_gf[0][i][2*j];

    // y-component
    real += pola_y_re(px,py,pz)*fnyquist_gf[1][i][2*j]   - pola_y_im(px,py,pz)*fnyquist_gf[1][i][2*j+1];
    imag += pola_y_re(px,py,pz)*fnyquist_gf[1][i][2*j+1] + pola_y_im(px,py,pz)*fnyquist_gf[1][i][2*j];

    // z-component
    real += pola_z_re(px,py,pz)*fnyquist_gf[2][i][2*j]   - pola_z_im(px,py,pz)*fnyquist_gf[2][i][2*j+1];
    imag += pola_z_re(px,py,pz)*fnyquist_gf[2][i][2*j+1] + pola_z_im(px,py,pz)*fnyquist_gf[2][i][2*j];

    // Add 3 directions
    amp = pw2(real) + pw2(imag);
  }

  return amp;
}


/* Function judging if do the cutoff process */
/*----------------------------------------------------------------------*/
bool cutoff_process_output(double pmagnitude, double px, double py, double pz)
{
  if(cutoff_method==0)
  {
    if(pmagnitude>modecutoff)
      return true;
    else
      return false;
  }
  else if(cutoff_method==1)
  {
    if(abs(px)>(double)N/4. || abs(py)>(double)N/4. || abs(pz)>(double)N/4.)
      return true;
    else
      return false;
  }
  else
  {
    printf("ERROR: Value of mutoff_method is wrong.\n");
  }
}

/* Power Spectra */
/*----------------------------------------------------------------------*/
// Inflaton(Axion) Spectra
void spectra_ax()
{
  // kでのスペクトル,k_effでのスペクトル,保存時刻
  static FILE *spectra_ax_,*spectra_eff_ax_,*spectratimes_;
  const int numbins=(int)(1.73205*(N/2))+1; // Number of bins (bin spacing=lattice spacing in Fourier space) = sqrt(NDIMS)*(N/2)+1. Set for 3D (i.e. biggest possible).
  int numpoints[numbins],numpoints_eff[numbins];
  double p[numbins],pave[numbins],peff_ave[numbins],f2[numbins];
  double p_eff[numbins],pave_eff[numbins],f2_eff[numbins];

  double pmagnitude; // pmagnitude = Sqrt(px^2+py^2+pz^2)
  double peffmagnitude;
  double dp=2.*pi/L; // Size of grid spacing in momentum space
  double fp2; // Square magnitude of field (fp2) and derivative (fpd2) for a given mode

  // parameter with dimension
  int i,j,k,px,py,pz; // px, py, and pz are components of momentum in units of grid spacing
  double fnyquist[N][2*N];
  int arraysize[]={N,N,N};
  double av_fd=0.;
  double xi=0.;


  static int first=1;
  if(first) // Open output files
  {
    sprintf(name_,"spectra_ax%s",ext_);
    spectra_ax_=fopen(name_,mode_);
    sprintf(name_,"spectra_eff_ax%s",ext_);
    spectra_eff_ax_=fopen(name_,mode_);
    sprintf(name_,"spectratimes%s",ext_);
    spectratimes_=fopen(name_,mode_);
    first=0;
  }

  // Calculate magnitude of momentum in each bin
  for(i=0;i<numbins;i++)
  {
    p[i]=dp*i;
    p_eff[i]=dp*i;
  }


    for(i=0;i<numbins;i++) // Initialize all bins to 0
    {
      numpoints[i]=0; // Number of points in the bin
      pave[i]=0.;
      peff_ave[i]=0.;
      f2[i]=0.; // |f_p|^2
      numpoints_eff[i]=0;
      pave_eff[i]=0.;
      f2_eff[i]=0.;
    }

    // FFT
    fftrn((double *)f[0],(double *)fnyquist,NDIMS,arraysize,1);


    // Loop runs over all gridpoints. All points with k<N/2 are in the array f, while points with k=N/2 are in the array fnyquist.
    // px and py go over all mode values in wrap-around order, rising from 0 to N/2 and then from -N/2+1 to -1
    for(i=0;i<N;i++)
    {
      px=(i<=N/2 ? i : i-N);
      for(j=0;j<N;j++)
      {
        py=(j<=N/2 ? j : j-N);
        // Modes with 0<k<N/2 are counted twice to account for the modes f(-p)=f(p)* that aren't explicitly included in the lattice
        for(k=1;k<N/2;k++) 
        {
          pz=k;
          fp2=pw2(f[0][i][j][2*k])+pw2(f[0][i][j][2*k+1]);

          // Lattice上の運動量
          pmagnitude=sqrt(pw2(px)+pw2(py)+pw2(pz)); // Magnitude of momentum of mode in units of momentum grid spacing
          peffmagnitude=(double)N*sqrt(pw2(sin(2.*pi*(double)px/(double)N))+pw2(sin(2.*pi*(double)py/(double)N))+pw2(sin(2.*pi*(double)pz/(double)N)))/(2.*pi);

          /* cutoff */
          if(cutoff_speax)
          {
            if(cutoff_process_output(pmagnitude,px,py,pz))
              continue;
          }

          numpoints[(int)pmagnitude] += 2;
          pave[(int)pmagnitude] += 2.*pmagnitude;
          peff_ave[(int)pmagnitude] += 2.*peffmagnitude;
          f2[(int)pmagnitude] += 2.*fp2;
    
          // 有効運動量
          numpoints_eff[(int)peffmagnitude] += 2;
          pave_eff[(int)peffmagnitude] += 2.*peffmagnitude;
          f2_eff[(int)peffmagnitude] += 2.*fp2;
        }

        // Modes with k=0 or k=N/2 are only counted once
        for(k=0;k<=N/2;k+=N/2)
        {
          pz=k;
          pmagnitude=sqrt(pw2(px)+pw2(py)+pw2(pz));

          /* cutoff */
          if(cutoff_speax)
          {
            if(cutoff_process_output(pmagnitude,px,py,pz))
              continue;
          }

          peffmagnitude=(double)N*sqrt(pw2(sin(2.*pi*(double)px/(double)N))+pw2(sin(2.*pi*(double)py/(double)N)))/(2.*pi);

          if(k==0) // Amplitudes of k=0 modes
          {   
            fp2=pw2(f[0][i][j][0])+pw2(f[0][i][j][1]);
          }
          else // Modes with k=N/2 are stored in the array fnyquist
          {
            fp2=pw2(fnyquist[i][2*j])+pw2(fnyquist[i][2*j+1]);
          }
          numpoints[(int)pmagnitude]++; // Iterate the count of points in this bin
          pave[(int)pmagnitude] += pmagnitude;
          peff_ave[(int)pmagnitude] += peffmagnitude;
          f2[(int)pmagnitude] += fp2; // Add the power of this mode to the bin


          if((i==0 || i==N/2) && (j==0 || j==N/2))
          {
              continue;
          }
          numpoints_eff[(int)peffmagnitude]++;
          pave_eff[(int)peffmagnitude] += peffmagnitude;
          f2_eff[(int)peffmagnitude] += fp2;
        }
      }
    }
    
    for(i=0;i<numbins;i++)
    {
      if(numpoints[i]>0) // Convert sums to averages. (numpoints[i] should always be greater than zero.)
      {
        pave[i] = pave[i]*dp/numpoints[i];
        peff_ave[i] = peff_ave[i]*dp/numpoints[i];
        f2[i] = f2[i]/numpoints[i];
      }
      if(numpoints_eff[i]>0) // Convert sums to averages. (numpoints[i] should always be greater than zero.)
      {
        pave_eff[i] = pave_eff[i]*dp/numpoints_eff[i];
        f2_eff[i] = f2_eff[i]/numpoints_eff[i];
      }
      // Output the momentum, number of points, omega, and calculated spectra for each bin
      fprintf(spectra_ax_,"%e %d %e %e %e\n",// %e %e %e\n",
        p[i],numpoints[i],pave[i],peff_ave[i],f2[i]);//,norm1*fd2[i],norm2*ndensity[i],norm2*rhodensity[i]);
      fprintf(spectra_eff_ax_,"%e %d %e %e\n",
        p_eff[i],numpoints_eff[i],pave_eff[i],f2_eff[i]);
    }
    fprintf(spectra_ax_,"\n");
    fflush(spectra_ax_);
    fprintf(spectra_eff_ax_,"\n");
    fflush(spectra_eff_ax_);

    // Inverse FFT
    fftrn((double *)f[0],(double *)fnyquist,NDIMS,arraysize,-1);


  LOOP
    av_fd += AXIOND;
  av_fd = av_fd/(double)gridsize;
  xi = g*av_fd*a/(2.*ad);
  fprintf(spectratimes_,"%f %f %e %f %f %f\n",t,a,ad,av_fd,xi,log(a));
  fflush(spectratimes_);


  //PBH mass algorithm
  /*
  if(spbh)
  {
    static int firstk = 1;
    if(firstk)
    {
      for(i=0;i<numbins;i++)
      {
        if(pave[i]>0.00001)
          paveList[i]=pave[i];
        if(peff_ave[i]>0.00001)
          peff_aveList[i]=peff_ave[i];
        if(pave_eff[i]>0.00001)
          pave_effList[i]=pave_eff[i];
      }
      firstk = 0;
    }
  }*/

  return;
}


// Gauge Spectra
void spectra_gf()
{
  static FILE *spectra_gf_,*spectra_eff_gf_;
  const int numbins=(int)(1.73205*(N/2))+1;
  int numpoints[numbins],numpoints_eff[numbins];
  double p[numbins],pave[numbins],peff_ave[numbins],f2posi[numbins],f2nega[numbins];
  double p_eff[numbins],pave_eff[numbins],f2posi_eff[numbins],f2nega_eff[numbins];

  double pmagnitude;
  double peffmagnitude;
  double dp=2.*pi/L;
  double fp2posi;
  double fp2nega;

  int i,j,k,px,py,pz;
  int arraysize[]={N,N,N};


  static int first=1;
  if(first)
  {
    sprintf(name_,"spectra_gf%s",ext_);
    spectra_gf_=fopen(name_,mode_);
    sprintf(name_,"spectra_eff_gf%s",ext_);
    spectra_eff_gf_=fopen(name_,mode_);
    first=0;
  }

  // Calculate magnitude of momentum in each bin
  for(i=0;i<numbins;i++)
  {
    p[i]=dp*i;
    p_eff[i]=dp*i;
  }
    
  for(i=0;i<numbins;i++) // Initialize all bins to 0
  {
    numpoints[i]=0;
    pave[i]=0.;
    peff_ave[i] = 0.;
    f2posi[i]=0.;
    f2nega[i] = 0.;
    numpoints_eff[i]=0;
    pave_eff[i]=0.;
    f2posi_eff[i]=0.;
    f2nega_eff[i] = 0.;
  }

  // Initialize fnyquist list
  for(i=0;i<N;i++)
    for(j=0;j<2*N;j++)
    {
      fnyquist_gf[0][i][j] = 0.;
      fnyquist_gf[1][i][j] = 0.;
      fnyquist_gf[2][i][j] = 0.;
    }

  // FFT A1~A3
  fftrn((double *)f[1],(double *)fnyquist_gf[0],NDIMS,arraysize,1);
  fftrn((double *)f[2],(double *)fnyquist_gf[1],NDIMS,arraysize,1);
  fftrn((double *)f[3],(double *)fnyquist_gf[2],NDIMS,arraysize,1);


    // Loop runs over all gridpoints. All points with k<N/2 are in the array f, while points with k=N/2 are in the array fnyquist.
    // px and py go over all mode values in wrap-around order, rising from 0 to N/2 and then from -N/2+1 to -1
    for(i=0;i<N;i++)
    {
      px=(i<=N/2 ? i : i-N);
      for(j=0;j<N;j++)
      {
        py=(j<=N/2 ? j : j-N);
        // Modes with 0<k<N/2 are counted twice to account for the modes f(-p)=f(p)* that aren't explicitly included in the lattice
        for(k=1;k<N/2;k++) 
        {
          pz=k;

          // Lattice上の運動量
          pmagnitude=sqrt(pw2(px)+pw2(py)+pw2(pz));

          /* cutoff */
          if(cutoff_spegf)
          {
            if(cutoff_process_output(pmagnitude,px,py,pz))
              continue;
          }

          peffmagnitude=(double)N*sqrt(pw2(sin(2.*pi*(double)px/(double)N))+pw2(sin(2.*pi*(double)py/(double)N))+pw2(sin(2.*pi*(double)pz/(double)N)))/(2.*pi);

          fp2posi=amplitude_gfposi(px,py,pz,i,j,k);
          fp2nega=amplitude_gfnega(px,py,pz,i,j,k);
          numpoints[(int)pmagnitude] += 2;
          pave[(int)pmagnitude] += 2.*pmagnitude;
          peff_ave[(int)pmagnitude] += 2.*peffmagnitude;
          f2posi[(int)pmagnitude] += 2.*fp2posi;
          f2nega[(int)pmagnitude] += 2.*fp2nega;
    
          // 有効運動量
          numpoints_eff[(int)peffmagnitude] += 2;
          pave_eff[(int)peffmagnitude] += 2.*peffmagnitude;
          f2posi_eff[(int)peffmagnitude] += 2.*fp2posi;
          f2nega_eff[(int)peffmagnitude] += 2.*fp2nega;
        }

        // Modes with k=0 or k=N/2 are only counted once
        for(k=0;k<=N/2;k+=N/2)
        {
          pz=k;
          pmagnitude=sqrt(pw2(px)+pw2(py)+pw2(pz));

          /* cutoff */
          if(cutoff_spegf)
          {
            if(cutoff_process_output(pmagnitude,px,py,pz))
              continue;
          }

          peffmagnitude=(double)N*sqrt(pw2(sin(2.*pi*(double)px/(double)N))+pw2(sin(2.*pi*(double)py/(double)N))+pw2(sin(2.*pi*(double)pz/(double)N)))/(2.*pi);

          if(k==0)
          {   
            fp2posi=amplitude_gfposi(px,py,0.,i,j,0);
            fp2nega=amplitude_gfnega(px,py,0.,i,j,0);
          }
          else
          {
            fp2posi=amplitude_gfposi_nyquist(px,py,i,j);
            fp2nega=amplitude_gfnega_nyquist(px,py,i,j);
          }
          numpoints[(int)pmagnitude]++;
          pave[(int)pmagnitude] += pmagnitude;
          peff_ave[(int)pmagnitude] += peffmagnitude;
          f2posi[(int)pmagnitude] += fp2posi;
          f2nega[(int)pmagnitude] += fp2nega;

          numpoints_eff[(int)peffmagnitude]++;
          pave_eff[(int)peffmagnitude] += peffmagnitude;
          f2posi_eff[(int)peffmagnitude] += fp2posi;
          f2nega_eff[(int)peffmagnitude] += fp2nega;
        }
      }
    }
    
    for(i=0;i<numbins;i++)
    {
      if(numpoints[i]>0) // Convert sums to averages. (numpoints[i] should always be greater than zero.)
      {
        pave[i] = pave[i]*dp/numpoints[i];
        peff_ave[i] = peff_ave[i]*dp/numpoints[i];
        f2posi[i] = f2posi[i]/numpoints[i];
        f2nega[i] = f2nega[i]/numpoints[i];
      }
      if(numpoints_eff[i]>0) // Convert sums to averages. (numpoints[i] should always be greater than zero.)
      {
        pave_eff[i] = pave_eff[i]*dp/numpoints_eff[i];
        f2posi_eff[i] = f2posi_eff[i]/numpoints_eff[i];
        f2nega_eff[i] = f2nega_eff[i]/numpoints_eff[i];
      }
      // Output the momentum, number of points, omega, and calculated spectra for each bin
      fprintf(spectra_gf_,"%e %d %e %e %e %e\n",// %e %e %e\n",
        p[i],numpoints[i],pave[i],peff_ave[i],f2posi[i],f2nega[i]);//,norm1*fd2[i],norm2*ndensity[i],norm2*rhodensity[i]);
      fprintf(spectra_eff_gf_,"%e %d %e %e %e\n",
        p_eff[i],numpoints_eff[i],pave_eff[i],f2posi_eff[i],f2nega_eff[i]);
    }
    fprintf(spectra_gf_,"\n");
    fflush(spectra_gf_);
    fprintf(spectra_eff_gf_,"\n");
    fflush(spectra_eff_gf_);

    // Transform field values back to position space.
    fftrn((double *)f[1],(double *)fnyquist_gf[0],NDIMS,arraysize,-1);
    fftrn((double *)f[2],(double *)fnyquist_gf[1],NDIMS,arraysize,-1);
    fftrn((double *)f[3],(double *)fnyquist_gf[2],NDIMS,arraysize,-1);

  return;
}



#endif