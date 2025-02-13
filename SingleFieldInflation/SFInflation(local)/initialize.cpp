/* This file implement important function "initialize()" */

#include "latticeeasy.hpp"
#include <random>


// Generating [0,1] random number using the Park-Miller algorithm
#define randa 16807
#define randm 2147483647
#define randq 127773
#define randr 2836
double rand_uniform(void)
{
  if(seed<1) return(0.33);
  static int i=0;
  static int next=seed;
  if(!(next>0))
  {
    printf("ERROR: Invalid seed used in random number function. Using seed=1\n");
    next=1;
  }
  if(i==0)
    for(i=1;i<100;i++)
      rand_uniform();
  next = randa*(next%randq) - randr*(next/randq);
  if(next<0) next += randm;
  return ((double)next/(double)randm);
}
#undef randa
#undef randm
#undef randq
#undef randr


// Generating Inflaton Mode Function
void set_mode(double keff, double *field, double *deriv, int real)
{
  double uni_X1,uni_Y1,uni_X2,uni_Y2;
  double gauss_X1,gauss_Y1,gauss_X2,gauss_Y2;
  double f1_re,f1_im,f2_re,f2_im;

  double f_norm;
  double amp,coeff1,coeff2;

  if(hankel) // HANKEL
  {
    double arg = keff/hubble_init;
    double nu = sqrt(2.25 - dvdf2()/pw2(hubble_init));

    double bessel_00 = std::cyl_bessel_j(nu,arg);
    double neuman_00 = std::cyl_neumann(nu,arg);
    double bessel_m1 = std::cyl_bessel_j(nu-1.,arg);
    double neuman_m1 = std::cyl_neumann(nu-1.,arg);

    amp = sqrt(pw2(bessel_00)+pw2(neuman_00));
    coeff1 = (nu-1.5)*hubble_init - keff*(bessel_00*bessel_m1+neuman_00*neuman_m1)/pw2(amp);
    coeff2 = keff*(bessel_00*neuman_m1-bessel_m1*neuman_00)/pw2(amp);

    f_norm  = m*pow(L,1.5)*sqrt(pi/hubble_init)/(2.*pow(dx,3));

  }
  else // EXP
  {
    f_norm  = m*pow(L,1.5)/(sqrt(2.)*pow(pw2(keff)+dvdf2(),1./4.)*pow(dx,3));
  }


  if(keff<0.001)
  {
    field[0] = 0.;
    field[1] = 0.;
    deriv[0] = 0.;
    deriv[1] = 0.;
  }
  else
  {
    uni_X1 = rand_uniform();
    uni_Y1 = rand_uniform();
    gauss_X1 = sqrt(log(1./uni_X1));
    gauss_Y1 = 2.*pi*uni_Y1;

    uni_X2 = rand_uniform();
    uni_Y2 = rand_uniform();
    gauss_X2 = sqrt(log(1./uni_X2));
    gauss_Y2 = 2.*pi*uni_Y2;

    if(hankel)
    {
      f1_re = cos(gauss_Y1)*gauss_X1*f_norm*amp;
      f1_im = sin(gauss_Y1)*gauss_X1*f_norm*amp;

      f2_re = cos(gauss_Y2)*gauss_X2*f_norm*amp;
      f2_im = sin(gauss_Y2)*gauss_X2*f_norm*amp;

      field[0] = (f1_re + f2_re)/sqrt(2.);
      field[1] = (f1_im + f2_im)/sqrt(2.);

      deriv[0] = coeff1*field[0] + coeff2*f_norm*amp*(sin(gauss_Y1)*gauss_X1 - sin(gauss_Y2)*gauss_X2)/sqrt(2.);
      deriv[1] = coeff1*field[1] - coeff2*f_norm*amp*(cos(gauss_Y1)*gauss_X1 - cos(gauss_Y2)*gauss_X2)/sqrt(2.);
    }
    else
    {
      f1_re = cos(gauss_Y1)*gauss_X1*f_norm;
      f1_im = sin(gauss_Y1)*gauss_X1*f_norm;

      f2_re = cos(gauss_Y2)*gauss_X2*f_norm;
      f2_im = sin(gauss_Y2)*gauss_X2*f_norm;

      field[0] = (f1_re + f2_re)/sqrt(2.);
      field[1] = (f1_im + f2_im)/sqrt(2.);

      deriv[0] = - ad*field[0] - f_norm*sqrt(pw2(keff)+1.)*(sin(gauss_Y1)*gauss_X1 - sin(gauss_Y2)*gauss_X2)/sqrt(2.);
      deriv[1] = - ad*field[1] + f_norm*sqrt(pw2(keff)+1.)*(cos(gauss_Y1)*gauss_X1 - cos(gauss_Y2)*gauss_X2)/sqrt(2.);
    }

    if(real==1)
    {
      field[1] = 0.;
      deriv[1] = 0.;
    }
  }

  return;
}



// Setting Initial Values
void initialize(ARG)
{
  int fld;
  double keff;
  double initial_field_values[nflds];
  double initial_field_derivs[nflds];
  FILE *old_grid_;

  int i,j,k,iconj,jconj;
  double px,py,pz;
  double fnyquist[N][2*N],fdnyquist[N][2*N];
  int arraysize[]={N,N,N};

  double (*f)[N][N][N] = flist[0];
  double (*fd)[N][N][N] = flist[1];

  // DEBUG: Courant condition dt/dx < 1/Sqrt(ndims)
  if(dt>dx/sqrt((double)NDIMS))
  {
    printf("ERROR: Time step too large. The ratio dt/dx is currently %f but for stability should never exceed 1/sqrt(%d) (%f)\n",dt/dx,NDIMS,1./sqrt((double)NDIMS));
    printf("Adjust dt to AT MOST %e, and preferably somewhat smaller than that\n",dx/sqrt((double)NDIMS));
    exit(1);
  }

  // Output initializations - Set values of nfldsout and ext_
  if(alt_extension[0]!='\0') // If an alternate extension was given use that instead of the default "_<run_number>.dat"
    sprintf(ext_,"%s",alt_extension);


  /* The case of using grid image */
  /*---------------------------------------------------------------------------------*/
  if(continue_run>0 && (old_grid_=fopen("grid.img","rb")))
  {
    printf("Previously generated grid image found. Reading in data...\n");
    fread(&run_number,sizeof(run_number),1,old_grid_);
    run_number++;
    fread(&t0,sizeof(t0),1,old_grid_);
    if(t0>=tf)
    {
      printf("A grid image file was found in this directory with values stored at t=%f. To continue that run set tf to a later time. To start a new run move or rename the file grid.img.\n",t0);
      exit(1);
    }
    fread(&a,sizeof(a),1,old_grid_);
    fread(&ad,sizeof(ad),1,old_grid_);
    fread(f,sizeof(double),nflds*gridsize,old_grid_);
    fread(fd,sizeof(double),nflds*gridsize,old_grid_);
    fclose(old_grid_);
    if(continue_run==1)
      sprintf(mode_,"a+");
    else if(alt_extension[0]=='\0')
      sprintf(ext_,"_%d.dat",run_number);
    printf("Data read. Resuming run at t=%f\n",t0);
    output_parameters();
    return;
  }


  /* The case of initializing without grid image */
  /*---------------------------------------------------------------------------------*/
  printf("OK: Generating initial conditions for new run at t=0\n");
  t0=0;
  run_number=0;

  for(fld=0;fld<nflds;fld++)
  {
    if(fld<(int)(sizeof initfield/sizeof(double)))
      initial_field_values[fld]=initfield[fld];
    else
      initial_field_values[fld]=0.;

    if(fld<(int)(sizeof initderivs/sizeof(double)))
      initial_field_derivs[fld]=initderivs[fld];
    else
      initial_field_derivs[fld]=0.;
  }

  // Set hubble parameter
  hubble_init = sqrt((pw2(initderivs[0])+2.*hubble_pote())/6.);
  if(!(hubble_init>=0.))
  {
    printf("ERROR: Bad initial Hubble constant\n");
    exit(1);
  }
  ad=hubble_init;

  // Set Inflaton Fluctuation
  for(fld=0;fld<nflds;fld++) 
  {
    for(i=0;i<N;i++)
    {
      px=(i<=N/2 ? i : i-N);
      iconj=(i==0 ? 0 : N-i);
      for(j=0;j<N;j++)
      {
        py=(j<=N/2 ? j : j-N);

        /* Mode with 0<k<N/2 */
        for(k=1;k<N/2;k++)
        {
          pz=k;
          keff=2.*sqrt(pw2(sin(pi*(double)px/(double)N))+pw2(sin(pi*(double)py/(double)N))+pw2(sin(pi*(double)pz/(double)N)))/dx;
          set_mode(keff,&f[fld][i][j][2*k],&fd[fld][i][j][2*k],0);
        }

        /* Mode with k=0,N/2 */
        if(j>N/2 || (i>N/2 && (j==0 || j==N/2)))
        {
          jconj=(j==0 ? 0 : N-j);
          // k=0
          keff=2.*sqrt(pw2(sin(pi*(double)px/(double)N))+pw2(sin(pi*(double)py/(double)N)))/dx;
          set_mode(keff,&f[fld][i][j][0],&fd[fld][i][j][0],0);
          f[fld][iconj][jconj][0]=f[fld][i][j][0];
          f[fld][iconj][jconj][1]=-f[fld][i][j][1];
          fd[fld][iconj][jconj][0]=fd[fld][i][j][0];
          fd[fld][iconj][jconj][1]=-fd[fld][i][j][1];
          // k=N/2
          keff=2.*sqrt(pw2(sin(pi*(double)px/(double)N))+pw2(sin(pi*(double)py/(double)N))+1.)/dx;
          set_mode(keff,&fnyquist[i][2*j],&fdnyquist[i][2*j],0);
          fnyquist[iconj][2*jconj]=fnyquist[i][2*j];
          fnyquist[iconj][2*jconj+1]=-fnyquist[i][2*j+1];
          fdnyquist[iconj][2*jconj]=fdnyquist[i][2*j];
          fdnyquist[iconj][2*jconj+1]=-fdnyquist[i][2*j+1];
        }
        else if((i==0 || i==N/2) && (j==0 || j==N/2)) // The 8 "corners" of the lattice are set to real values
        {
          // k=0
          keff=2.*sqrt(pw2(sin(pi*(double)px/(double)N))+pw2(sin(pi*(double)py/(double)N)))/dx;
          if(keff>0.)
            set_mode(keff,&f[fld][i][j][0],&fd[fld][i][j][0],1);
          // k=N/2
          keff=2.*sqrt(pw2(sin(pi*(double)px/(double)N))+pw2(sin(pi*(double)py/(double)N))+1.)/dx;
          set_mode(keff,&fnyquist[i][2*j],&fdnyquist[i][2*j],1);
        }
      }
    }

    // Set zero mode
    f[fld][0][0][0]=0.;
    f[fld][0][0][1]=0.;
    fd[fld][0][0][0]=0.;
    fd[fld][0][0][1]=0.;

    // DEBUG: zero mode
    if(f[0][0][0][0]==0. && f[0][0][0][1]==0. && fd[0][0][0][0]==0. && fd[0][0][0][1]==0.)
      printf("OK: Zero mode is zero\n");
    else
      printf("ERROR: Zero mode is not zero\n");

    // Inverse FFT
    if(seed>=0)
    {
      fftrn((double *)f[fld],(double *)fnyquist,NDIMS,arraysize,-1);
      fftrn((double *)fd[fld],(double *)fdnyquist,NDIMS,arraysize,-1);
    }

    // DEBUG: if fluctuation is small
    if(f[0][0][0][0]<0.01 && f[0][1][50][32]<0.01)
      printf("OK: Fluctuation is small\n");
    else
      printf("ERROR: Flucuation is not small\n");

    // Adding back ground
    LOOP
    {
      FIELD(fld) += initial_field_values[fld];
      FIELDD(fld) += initial_field_derivs[fld];
    }
  }

  // Post processing
  save(0,flist);
  output_parameters();
  printf("OK: Finished initial conditions\n");

  return;
}







