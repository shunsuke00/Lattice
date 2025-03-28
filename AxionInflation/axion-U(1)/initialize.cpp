/* This file implements important function "initialize_ax()" and "initialize_gf()" */

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


// Generating Inflaton(Axion) Mode Function
inline void set_mode_ax(double keff, double *field, double *deriv, int real, bool cutoff)
{
  double uni_X,uni_Y;
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

    amp    = sqrt(pw2(bessel_00)+pw2(neuman_00));
    coeff1 = (nu-1.5)*hubble_init - keff*(bessel_00*bessel_m1+neuman_00*neuman_m1)/pw2(amp);
    coeff2 = keff*(bessel_00*neuman_m1-bessel_m1*neuman_00)/pw2(amp);

    f_norm  = m*pow(L,1.5)*sqrt(pi/hubble_init)/(2.*pow(dx,3));
  }
  else // EXP
  {
    f_norm = m*pow(L,1.5)/(sqrt(2.)*pow(pw2(keff)+dvdf2(),1./4.)*pow(dx,3.));
  }


  if(cutoff)
  {
    field[0] = 0.;
    field[1] = 0.;
    deriv[0] = 0.;
    deriv[1] = 0.;
  }
  else if(keff<0.01)
  {
    field[0] = 0.;
    field[1] = 0.;
    deriv[0] = 0.;
    deriv[1] = 0.;
  }
  else
  {
    uni_X = rand_uniform();
    uni_Y = rand_uniform();
    gauss_X1 = sqrt(log(1./uni_X));
    gauss_Y1 = 2.*pi*uni_Y;

    uni_X = rand_uniform();
    uni_Y = rand_uniform();
    gauss_X2 = sqrt(log(1./uni_X));
    gauss_Y2 = 2.*pi*uni_Y;

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

      deriv[0] = - ad*field[0] - f_norm*sqrt(pw2(keff)+dvdf2())*(sin(gauss_Y1)*gauss_X1 - sin(gauss_Y2)*gauss_X2)/sqrt(2.);
      deriv[1] = - ad*field[1] + f_norm*sqrt(pw2(keff)+dvdf2())*(cos(gauss_Y1)*gauss_X1 - cos(gauss_Y2)*gauss_X2)/sqrt(2.);
    }

    if(real==1)
    {
      field[1] = 0.;
      deriv[1] = 0.;
    }

  }

  return;
}



// Generating Gauge Mode Function [A_k]_mu
inline void set_mode_gf(double px, double py, double pz, double keff, double grow_amp, double grow_fac, double stab_amp, double stab_fac, int i, int j, int k, int real, bool cutoff)
{
  double uniX_po,uniY_po,uniX_ne,uniY_ne;
  double gaussX_po1,gaussY_po1,gaussX_ne1,gaussY_ne1;
  double gaussX_po2,gaussY_po2,gaussX_ne2,gaussY_ne2;

  double u_po_re,u_po_im,u_ne_re,u_ne_im;
  double du_po_re,du_po_im,du_ne_re,du_ne_im;

  double *field;
  double *deriv;
  int n;

  double xi = g*initderiv[0]/(2.*hubble_init);
  double f_norm = m*pow(L,1.5)/(pow(dx,3.)*sqrt(2.*keff));


  uniX_po = rand_uniform();
  uniY_po = rand_uniform();
  uniX_ne = rand_uniform();
  uniY_ne = rand_uniform();
  gaussX_po1 = sqrt(log(1./uniX_po));
  gaussY_po1 = 2.*pi*uniY_po;
  gaussX_ne1 = sqrt(log(1./uniX_ne));
  gaussY_ne1 = 2.*pi*uniY_ne;

  uniX_po = rand_uniform();
  uniY_po = rand_uniform();
  uniX_ne = rand_uniform();
  uniY_ne = rand_uniform();
  gaussX_po2 = sqrt(log(1./uniX_po));
  gaussY_po2 = 2.*pi*uniY_po;
  gaussX_ne2 = sqrt(log(1./uniX_ne));
  gaussY_ne2 = 2.*pi*uniY_ne;


  // + Polarization 
  if(coulomb_po) // Coulomb Wave Function
  {
    u_po_re = f_norm*stab_amp*(gaussX_po1*cos(gaussY_po1)+gaussX_po2*cos(gaussY_po2))/sqrt(2.);
    u_po_im = f_norm*stab_amp*(gaussX_po1*sin(gaussY_po1)+gaussX_po2*sin(gaussY_po2))/sqrt(2.);

    du_po_re = -(hubble_init-keff*xi-keff*sqrt(1.+pw2(xi))*stab_fac)*u_po_re + keff*f_norm*(gaussX_po1*sin(gaussY_po1)-gaussX_po2*sin(gaussY_po2))/(stab_amp*sqrt(2.));
    du_po_im = -(hubble_init-keff*xi-keff*sqrt(1.+pw2(xi))*stab_fac)*u_po_im - keff*f_norm*(gaussX_po1*cos(gaussY_po1)-gaussX_po2*cos(gaussY_po2))/(stab_amp*sqrt(2.));
  }
  else  // EXP
  {
    u_po_re = f_norm*(gaussX_po1*cos(gaussY_po1)+gaussX_po2*cos(gaussY_po2))/sqrt(2.);
    u_po_im = f_norm*(gaussX_po1*sin(gaussY_po1)+gaussX_po2*sin(gaussY_po2))/sqrt(2.);

    du_po_re = f_norm*keff*(-gaussX_po1*sin(gaussY_po1)+gaussX_po2*sin(gaussY_po2))/sqrt(2.);
    du_po_im = f_norm*keff*(+gaussX_po1*cos(gaussY_po1)-gaussX_po2*cos(gaussY_po2))/sqrt(2.);
  }

  // - Polarization
  if(coulomb_ne) // Coulomb Wave Function
  {
    u_ne_re = f_norm*grow_amp*(gaussX_ne1*cos(gaussY_ne1)+gaussX_ne2*cos(gaussY_ne2))/sqrt(2.);
    u_ne_im = f_norm*grow_amp*(gaussX_ne1*sin(gaussY_ne1)+gaussX_ne2*sin(gaussY_ne2))/sqrt(2.);

    du_ne_re = -(hubble_init+keff*xi-keff*sqrt(1.+pw2(xi))*grow_fac)*u_ne_re + keff*f_norm*(gaussX_ne1*sin(gaussY_ne1)-gaussX_ne2*sin(gaussY_ne2))/(grow_amp*sqrt(2.));
    du_ne_im = -(hubble_init+keff*xi-keff*sqrt(1.+pw2(xi))*grow_fac)*u_ne_im - keff*f_norm*(gaussX_ne1*cos(gaussY_ne1)-gaussX_ne2*cos(gaussY_ne2))/(grow_amp*sqrt(2.));
  }
  else  // EXP
  {
    u_ne_re = f_norm*(gaussX_ne1*cos(gaussY_ne1)+gaussX_ne2*cos(gaussY_ne2))/sqrt(2.);
    u_ne_im = f_norm*(gaussX_ne1*sin(gaussY_ne1)+gaussX_ne2*sin(gaussY_ne2))/sqrt(2.);

    du_ne_re = f_norm*keff*(-gaussX_ne1*sin(gaussY_ne1)+gaussX_ne2*sin(gaussY_ne2))/sqrt(2.);
    du_ne_im = f_norm*keff*(+gaussX_ne1*cos(gaussY_ne1)-gaussX_ne2*cos(gaussY_ne2))/sqrt(2.);
  }


  if(cutoff)
  {
    for(n=2;n<=4;n++)
    {
      field = &f[n][i][j][2*k];
      deriv = &fd[n][i][j][2*k];
      field[0] = 0.;
      field[1] = 0.;
      deriv[0] = 0.;
      deriv[1] = 0.;
    }
  }
  else if(keff < 0.01)//zero mode
  {
    for(n=2;n<=4;n++)
    {
      field = &f[n][i][j][2*k];
      deriv = &fd[n][i][j][2*k];
      field[0] = 0.;
      field[1] = 0.;
      deriv[0] = 0.;
      deriv[1] = 0.;
    }
  }
  else
  {
    if(px==0. && py==0.)
    {
      // x-component
      field = &f[2][i][j][2*k];
      deriv = &fd[2][i][j][2*k];

      field[0] = (u_po_re+u_ne_re)/sqrt(2.);
      field[1] = (u_po_im+u_ne_im)/sqrt(2.);
      deriv[0] = (du_po_re+du_ne_re)/sqrt(2.);
      deriv[1] = (du_po_im+du_ne_im)/sqrt(2.);

      // y-component
      field = &f[3][i][j][2*k];
      deriv = &fd[3][i][j][2*k];

      field[0] = (-u_po_im+u_ne_im)/sqrt(2.);
      field[1] = (u_po_re-u_ne_re)/sqrt(2.);
      deriv[0] = (-du_po_re+du_ne_re)/sqrt(2.);
      deriv[1] = (du_po_re-du_ne_re)/sqrt(2.);

      // z-component
      field = &f[4][i][j][2*k];
      deriv = &fd[4][i][j][2*k];

      field[0] = 0.;
      field[1] = 0.;
      deriv[0] = 0.;
      deriv[1] = 0.;
    }
    else
    {
      // x-component
      field = &f[2][i][j][2*k];
      deriv = &fd[2][i][j][2*k];

      field[0] = pola_x_re(px,py,pz)*(u_po_re+u_ne_re) - pola_x_im(px,py,pz)*(u_po_im-u_ne_im);
      field[1] = pola_x_re(px,py,pz)*(u_po_im+u_ne_im) + pola_x_im(px,py,pz)*(u_po_re-u_ne_re);
      deriv[0] = pola_x_re(px,py,pz)*(du_po_re+du_ne_re) - pola_x_im(px,py,pz)*(du_po_im-du_ne_im);
      deriv[1] = pola_x_re(px,py,pz)*(du_po_im+du_ne_im) + pola_x_im(px,py,pz)*(du_po_re-du_ne_re);

      // y-component
      field = &f[3][i][j][2*k];
      deriv = &fd[3][i][j][2*k];

      field[0] = pola_y_re(px,py,pz)*(u_po_re+u_ne_re) - pola_y_im(px,py,pz)*(u_po_im-u_ne_im);
      field[1] = pola_y_re(px,py,pz)*(u_po_im+u_ne_im) + pola_y_im(px,py,pz)*(u_po_re-u_ne_re);
      deriv[0] = pola_y_re(px,py,pz)*(du_po_re+du_ne_re) - pola_y_im(px,py,pz)*(du_po_im-du_ne_im);
      deriv[1] = pola_y_re(px,py,pz)*(du_po_im+du_ne_im) + pola_y_im(px,py,pz)*(du_po_re-du_ne_re);

      // z-component
      field = &f[4][i][j][2*k];
      deriv = &fd[4][i][j][2*k];

      field[0] = pola_z_re(px,py,pz)*(u_po_re+u_ne_re) - pola_z_im(px,py,pz)*(u_po_im-u_ne_im);
      field[1] = pola_z_re(px,py,pz)*(u_po_im+u_ne_im) + pola_z_im(px,py,pz)*(u_po_re-u_ne_re);
      deriv[0] = pola_z_re(px,py,pz)*(du_po_re+du_ne_re) - pola_z_im(px,py,pz)*(du_po_im-du_ne_im);
      deriv[1] = pola_z_re(px,py,pz)*(du_po_im+du_ne_im) + pola_z_im(px,py,pz)*(du_po_re-du_ne_re);
    }
  }

  if(real==1)
  {
    for(n=2;n<=4;n++)
    {
      field = &f[n][i][j][2*k];
      deriv = &fd[n][i][j][2*k];
      field[1] = 0.;
      deriv[1] = 0.;
    }
  }
}


inline void set_mode_gf_nyquist(double px, double py, double keff, double grow_amp, double grow_fac, double stab_amp, double stab_fac, int i, int j, int real, bool cutoff)
{
  double pz = (double)N/2.;
  int k = N/2;

  double uniX_po,uniY_po,uniX_ne,uniY_ne;
  double gaussX_po1,gaussY_po1,gaussX_ne1,gaussY_ne1;
  double gaussX_po2,gaussY_po2,gaussX_ne2,gaussY_ne2;
  
  double u_po_re,u_po_im,u_ne_re,u_ne_im;
  double du_po_re,du_po_im,du_ne_re,du_ne_im;

  double *field;
  double *deriv;
  int n;

  double xi = g*initderiv[0]/(2.*hubble_init);
  double f_norm = m*pow(L,1.5)/(pow(dx,3.)*sqrt(2.*keff));


  uniX_po = rand_uniform();
  uniY_po = rand_uniform();
  uniX_ne = rand_uniform();
  uniY_ne = rand_uniform();
  gaussX_po1 = sqrt(log(1./uniX_po));
  gaussY_po1 = 2.*pi*uniY_po;
  gaussX_ne1 = sqrt(log(1./uniX_ne));
  gaussY_ne1 = 2.*pi*uniY_ne;

  uniX_po = rand_uniform();
  uniY_po = rand_uniform();
  uniX_ne = rand_uniform();
  uniY_ne = rand_uniform();
  gaussX_po2 = sqrt(log(1./uniX_po));
  gaussY_po2 = 2.*pi*uniY_po;
  gaussX_ne2 = sqrt(log(1./uniX_ne));
  gaussY_ne2 = 2.*pi*uniY_ne;


  // + Polarization 
  if(coulomb_po) // Coulomb Wave Function
  {
    u_po_re = f_norm*stab_amp*(gaussX_po1*cos(gaussY_po1)+gaussX_po2*cos(gaussY_po2))/sqrt(2.);
    u_po_im = f_norm*stab_amp*(gaussX_po1*sin(gaussY_po1)+gaussX_po2*sin(gaussY_po2))/sqrt(2.);

    du_po_re = -(hubble_init-keff*xi-keff*sqrt(1.+pw2(xi))*stab_fac)*u_po_re + keff*f_norm*(gaussX_po1*sin(gaussY_po1)-gaussX_po2*sin(gaussY_po2))/(stab_amp*sqrt(2.));
    du_po_im = -(hubble_init-keff*xi-keff*sqrt(1.+pw2(xi))*stab_fac)*u_po_im - keff*f_norm*(gaussX_po1*cos(gaussY_po1)-gaussX_po2*cos(gaussY_po2))/(stab_amp*sqrt(2.));
  }
  else  // EXP
  {
    u_po_re = f_norm*(gaussX_po1*cos(gaussY_po1)+gaussX_po2*cos(gaussY_po2))/sqrt(2.);
    u_po_im = f_norm*(gaussX_po1*sin(gaussY_po1)+gaussX_po2*sin(gaussY_po2))/sqrt(2.);

    du_po_re = f_norm*keff*(-gaussX_po1*sin(gaussY_po1)+gaussX_po2*sin(gaussY_po2))/sqrt(2.);
    du_po_im = f_norm*keff*(+gaussX_po1*cos(gaussY_po1)-gaussX_po2*cos(gaussY_po2))/sqrt(2.);
  }

  // - Polarization
  if(coulomb_ne) // Coulomb Wave Function
  {
    u_ne_re = f_norm*grow_amp*(gaussX_ne1*cos(gaussY_ne1)+gaussX_ne2*cos(gaussY_ne2))/sqrt(2.);
    u_ne_im = f_norm*grow_amp*(gaussX_ne1*sin(gaussY_ne1)+gaussX_ne2*sin(gaussY_ne2))/sqrt(2.);

    du_ne_re = -(hubble_init+keff*xi-keff*sqrt(1.+pw2(xi))*grow_fac)*u_ne_re + keff*f_norm*(gaussX_ne1*sin(gaussY_ne1)-gaussX_ne2*sin(gaussY_ne2))/(grow_amp*sqrt(2.));
    du_ne_im = -(hubble_init+keff*xi-keff*sqrt(1.+pw2(xi))*grow_fac)*u_ne_im - keff*f_norm*(gaussX_ne1*cos(gaussY_ne1)-gaussX_ne2*cos(gaussY_ne2))/(grow_amp*sqrt(2.));
  }
  else  // EXP
  {
    u_ne_re = f_norm*(gaussX_ne1*cos(gaussY_ne1)+gaussX_ne2*cos(gaussY_ne2))/sqrt(2.);
    u_ne_im = f_norm*(gaussX_ne1*sin(gaussY_ne1)+gaussX_ne2*sin(gaussY_ne2))/sqrt(2.);

    du_ne_re = f_norm*keff*(-gaussX_ne1*sin(gaussY_ne1)+gaussX_ne2*sin(gaussY_ne2))/sqrt(2.);
    du_ne_im = f_norm*keff*(+gaussX_ne1*cos(gaussY_ne1)-gaussX_ne2*cos(gaussY_ne2))/sqrt(2.);
  }


  if(cutoff)
  {
    for(n=0;n<=2;n++)
    {
      field = &fnyquist_gf[n][i][2*j];
      deriv = &fdnyquist_gf[n][i][2*j];
      field[0] = 0.;
      field[1] = 0.;
      deriv[0] = 0.;
      deriv[1] = 0.;
    }
  }
  else if(keff < 0.01)//zero mode
  {
    for(n=0;n<=2;n++)
    {
      field = &fnyquist_gf[n][i][2*j];
      deriv = &fdnyquist_gf[n][i][2*j];
      field[0] = 0.;
      field[1] = 0.;
      deriv[0] = 0.;
      deriv[1] = 0.;
    }
  }
  else//general mode
  {
    // x-component
    field = &fnyquist_gf[0][i][2*j];
    deriv = &fdnyquist_gf[0][i][2*j];

    field[0] = pola_x_re(px,py,pz)*(u_po_re+u_ne_re) - pola_x_im(px,py,pz)*(u_po_im-u_ne_im);
    field[1] = pola_x_re(px,py,pz)*(u_po_im+u_ne_im) + pola_x_im(px,py,pz)*(u_po_re-u_ne_re);
    deriv[0] = pola_x_re(px,py,pz)*(du_po_re+du_ne_re) - pola_x_im(px,py,pz)*(du_po_im-du_ne_im);
    deriv[1] = pola_x_re(px,py,pz)*(du_po_im+du_ne_im) + pola_x_im(px,py,pz)*(du_po_re-du_ne_re);

    // y-component
    field = &fnyquist_gf[1][i][2*j];
    deriv = &fdnyquist_gf[1][i][2*j];

    field[0] = pola_y_re(px,py,pz)*(u_po_re+u_ne_re) - pola_y_im(px,py,pz)*(u_po_im-u_ne_im);
    field[1] = pola_y_re(px,py,pz)*(u_po_im+u_ne_im) + pola_y_im(px,py,pz)*(u_po_re-u_ne_re);
    deriv[0] = pola_y_re(px,py,pz)*(du_po_re+du_ne_re) - pola_y_im(px,py,pz)*(du_po_im-du_ne_im);
    deriv[1] = pola_y_re(px,py,pz)*(du_po_im+du_ne_im) + pola_y_im(px,py,pz)*(du_po_re-du_ne_re);

    // z-component
    field = &fnyquist_gf[2][i][2*j];
    deriv = &fdnyquist_gf[2][i][2*j];
    
    field[0] = pola_z_re(px,py,pz)*(u_po_re+u_ne_re) - pola_z_im(px,py,pz)*(u_po_im-u_ne_im);
    field[1] = pola_z_re(px,py,pz)*(u_po_im+u_ne_im) + pola_z_im(px,py,pz)*(u_po_re-u_ne_re);
    deriv[0] = pola_z_re(px,py,pz)*(du_po_re+du_ne_re) - pola_z_im(px,py,pz)*(du_po_im-du_ne_im);
    deriv[1] = pola_z_re(px,py,pz)*(du_po_im+du_ne_im) + pola_z_im(px,py,pz)*(du_po_re-du_ne_re);
  }

  if(real==1)
  {
    for(n=0;n<=2;n++)
    {
      field = &fnyquist_gf[n][i][2*j];
      deriv = &fdnyquist_gf[n][i][2*j];
      field[1] = 0.;
      deriv[1] = 0.;
    }
  }
}


// Cutoff process
bool cutoff_process_initial(double px, double py, double pz)
{
  if(cutoff_method==0)
    return (sqrt(pw2(px) + pw2(py) + pw2(pz)) > modecutoff);
  else if(cutoff_method==1)
    return (abs(px)>(double)N/4. || abs(py)>(double)N/4. || abs(pz)>(double)N/4.);
  else
  {
    printf("ERROR: Value of mutoff_method is wrong.\n");
    return true;
  }
}


// Set initial parameters and axion
void initialize_ax()
{
  double keff;
  FILE *old_grid_; // Used for reading in previously generated data as initial conditions

  int i,j,k,iconj,jconj; // iconj and jconj are the indices where conjugate modes to i and j are stored
  double px,py,pz; // Components of momentum in units of grid spacing
  double fnyquist[N][2*N],fdnyquist[N][2*N];
  int arraysize[]={N,N,N};
  bool cutoff;

  // DEBUG: Courant condition, dt/dx < 1/Sqrt(ndims)
  if(dt>dx/sqrt((double)NDIMS))
  {
    printf("ERROR: Time step too large. The ratio dt/dx is currently %f but for stability should never exceed 1/sqrt(%d) (%f)\n",dt/dx,NDIMS,1./sqrt((double)NDIMS));
    printf("Adjust dt to AT MOST %e, and preferably somewhat smaller than that.\n",dx/sqrt((double)NDIMS));
    exit(1);
  }

  // Output initializations - Set values of nfldsout and ext_
  if(alt_extension[0]!='\0') // If an alternate extension was given use that instead of the default "_<run_number>.dat"
    sprintf(ext_,"%s",alt_extension);


  /* The case of using grid image */
  /*----------------------------------------------------------------*/
  if(continue_run>0 && (old_grid_=fopen("grid.img","rb")))
  {
    printf("Previously generated grid image found. Reading in data...\n");
    fread(&run_number,sizeof(run_number),1,old_grid_);
    run_number++;
    fread(&t0,sizeof(t0),1,old_grid_);
    if(t0>=tf) // Check to make sure that the time in the old grid is earlier than the final time for this run
    {
      printf("A grid image file was found in this directory with values stored at t=%f. To continue that run set tf to a later time. To start a new run move or rename the file grid.img.\n",t0);
      exit(1);
    }
    fread(&a,sizeof(a),1,old_grid_);
    fread(&ad,sizeof(ad),1,old_grid_);
    fread(f,sizeof(double),nflds*gridsize,old_grid_);
    fread(fd,sizeof(double),nflds*gridsize,old_grid_);
    fclose(old_grid_);
    if(continue_run==1) // Option to append new data to old data files
      sprintf(mode_,"a+");
    else if(alt_extension[0]=='\0') // If no alternate extension was given set filename extension to indicate run number
      sprintf(ext_,"_%d.dat",run_number);
    printf("Data read. Resuming run at t=%f\n",t0);
    output_parameters();
    return;
  }


  /* The case of initializing without grid image */
  /*----------------------------------------------------------------*/
  printf("OK: Generating initial conditions for new run at t=0\n");
  t0=0.;
  run_number=0;

  // Set Hubble parameter
  hubble_init = sqrt((pw2(initderiv[0])+2.*hubble_pote())/6.);
  if(!(hubble_init>=0.))
  {
    printf("ERROR: Bad initial Hubble constant. Exiting.\n");
    exit(1);
  }
  ad=hubble_init;


  // Set Axion fluctuation
  // px,py,pz = from 0 to N/2 and then from -N/2+1 to -1
  for(i=0;i<N;i++)
  {
    px=(i<=N/2 ? i : i-N);
    iconj=(i==0 ? 0 : N-i);
    for(j=0;j<N;j++)
    {
      py=(j<=N/2 ? j : j-N);

      // Set modes for 0<k<N/2
      for(k=1;k<N/2;k++)
      {
        pz=k;

        if(cutoff_iniax)
          cutoff = cutoff_process_initial(px,py,pz);
        else
          cutoff = false;

        keff=sqrt(pw2(sin(2.*pi*(double)px/(double)N))+pw2(sin(2.*pi*(double)py/(double)N))+pw2(sin(2.*pi*(double)pz/(double)N)))/dx;
        set_mode_ax(keff,&f[0][i][j][2*k],&fd[0][i][j][2*k],0,cutoff);
      }

      // Set modes with k=0 or N/2.
      if(j>N/2 || (i>N/2 && (j==0 || j==N/2)))
      {
        jconj=(j==0 ? 0 : N-j);

        // k=0
        pz = 0.;
        if(cutoff_iniax)
          cutoff = cutoff_process_initial(px,py,pz);
        else
          cutoff = false;

        keff=sqrt(pw2(sin(2.*pi*(double)px/(double)N))+pw2(sin(2.*pi*(double)py/(double)N)))/dx;
        set_mode_ax(keff,&f[0][i][j][0],&fd[0][i][j][0],0,cutoff);
        f[0][iconj][jconj][0]=f[0][i][j][0];
        f[0][iconj][jconj][1]=-f[0][i][j][1];
        fd[0][iconj][jconj][0]=fd[0][i][j][0];
        fd[0][iconj][jconj][1]=-fd[0][i][j][1];

        // k=N/2
        pz = (double)N/2.;
        if(cutoff_iniax)
          cutoff = cutoff_process_initial(px,py,pz);
        else
          cutoff = false;

        keff=sqrt(pw2(sin(2.*pi*(double)px/(double)N))+pw2(sin(2.*pi*(double)py/(double)N)))/dx;
        set_mode_ax(keff,&fnyquist[i][2*j],&fdnyquist[i][2*j],0,cutoff);
        fnyquist[iconj][2*jconj]=fnyquist[i][2*j];
        fnyquist[iconj][2*jconj+1]=-fnyquist[i][2*j+1];
        fdnyquist[iconj][2*jconj]=fdnyquist[i][2*j];
        fdnyquist[iconj][2*jconj+1]=-fdnyquist[i][2*j+1];
      }
      else if((i==0 || i==N/2) && (j==0 || j==N/2)) // The 8 "corners" of the lattice are set to real values
      {
        cutoff = false;
        keff=0.;
        if(i==0 && j==0)
        {
          set_mode_ax(keff,&f[0][i][j][0],&fd[0][i][j][0],0,cutoff);
        }
        else
        {
          set_mode_ax(keff,&f[0][i][j][0],&fd[0][i][j][0],1,cutoff);
        }
        keff=0.;
        set_mode_ax(keff,&fnyquist[i][2*j],&fdnyquist[i][2*j],1,cutoff);
      }
    }
  }

  // DEBUG: zero mode
  if(f[0][0][0][0]==0. && f[0][0][0][1]==0. && fd[0][0][0][0]==0. && fd[0][0][0][1]==0.)
    printf("OK: Axion zero mode is 0 (i=j=k=0 mode)\n");
  else
  {
    printf("ERROR: Axion zero mode is large\n");
    exit(1);
  }

  // Inverse FFT
  if(seed>=0)
  {
    fftrn((double *)f[0],(double *)fnyquist,NDIMS,arraysize,-1);
    fftrn((double *)fd[0],(double *)fdnyquist,NDIMS,arraysize,-1);
  }

  // DEBUG: whether fluctuation is small
  if(f[0][0][0][0]<0.01 && f[0][1][50][49]<0.01)
    printf("OK: Axion fluctuation is small\n");
  else
  {
    printf("ERROR: Axion fluctuation is large\n");
    exit(1);
  }

  // Add back ground
  LOOP
  {
    AXION  += initfield[0];
    AXIOND += initderiv[0];
  }

  // Post processing
  printf("OK: Finished initial conditions of Inflaton\n");
  return;
}


// Set gauge field
void initialize_gf()
{
  double keff;
  int fld;

  std::fstream growAmp,growFac,stabAmp,stabFac;
  std::string path_cwf = path_base;
  std::string path_size;
#if BACKREACTION==0
  path_cwf += "weak/";
  if(L==2.0)
    path_size = "2/";
  else if(L==4.0)
    path_size = "4/";
  else if(L==1.0)
    path_size = "1/";
  else
  {
    printf("ERROR: Coulomb file name is wrong\n");
    exit(1);
  }
#elif BACKREACTION==1
  path_cwf += "strong/";
  if(L==1.5)
    path_size = "1.5/";
  else if(L==0.4)
    path_size = "0.4/";
  else if(L==0.2)
    path_size = "0.2/";
  else
  {
    printf("ERROR: Coulomb file name is wrong\n");
    exit(1);
  }
#elif BACKREACTION==2
  path_cwf += "bumpy/";
  if(L==0.5)
    path_size = "0.5/";
  else if(L==1.2)
  {
    if(initfield[0]==-4.5042820)
      path_size = "H=1.93/L=1.2/";
    else if(initfield[0]==-4.5929167)
      path_size = "H=1.939/L=1.2/";
    else
      path_size = "1.2/";
  }
  else if(L==3.0)
  {
    if(initfield[0]==-4.6205366)
      path_size = "H=1.941/L=3/";
    else if(initfield[0]==-4.5929167)
      path_size = "H=1.939/L=3/";
    else if(initfield[0]==-4.5042820)
      path_size = "H=1.93/L=3/";
    else if(initfield[0]==-4.3387189)
      path_size = "H=1.9/L=3/";
  }
  else if(L==3.2)
  {
    if(initfield[0]==-4.5042820)
      path_size = "H=1.93/L=3.2/";
    else if(initfield[0]==-4.4273861)
      path_size = "H=1.92/L=3.2/";
  }
  else if(L==3.5)
  {
    if(initfield[0]==-4.4273861)
      path_size = "H=1.92/L=3.5/";
  }
  else
  {
    printf("ERROR: Coulomb file name is wrong\n");
    exit(1);
  }
#endif
  path_cwf += std::to_string(N) + "/" + path_size;

  std::string path_growAmp = path_cwf + "growAmp.dat";
  std::string path_growFac = path_cwf + "growFac.dat";
  std::string path_stabAmp = path_cwf + "stabAmp.dat";
  std::string path_stabFac = path_cwf + "stabFac.dat";

  int i,j,k,iconj,jconj; // iconj and jconj are the indices where conjugate modes to i and j are stored
  double px,py,pz;
  int arraysize[]={N,N,N};
  bool cutoff;
  double grow_amp,grow_fac,stab_amp,stab_fac;


  /* Gauge fluctuation */
  /*---------------------------------------------------------*/
  // File open
  growAmp.open(path_growAmp,std::ios::in);
  growFac.open(path_growFac,std::ios::in);
  stabAmp.open(path_stabAmp,std::ios::in);
  stabFac.open(path_stabFac,std::ios::in);
  
  if(!growAmp.is_open() || !growFac.is_open() || !stabAmp.is_open() || !stabFac.is_open())
  {
    printf("ERROR: coulomb files are not opened.\n");
    exit(1);
  }
  else
  {
    printf("OK: coulomb files are opened.\n");
    std::cout << "path1: " << path_growAmp.substr(50) << std::endl;
    std::cout << "path2: " << path_growFac.substr(50) << std::endl;
    std::cout << "path3: " << path_stabAmp.substr(50) << std::endl;
    std::cout << "path4: " << path_stabFac.substr(50) << std::endl;
  }

  // Substitude cwf values with for loop
  for(i=0;i<N;i++)
  {
    px=(i<=N/2 ? i : i-N);
    iconj=(i==0 ? 0 : N-i);
    for(j=0;j<N;j++)
    {
      py=(j<=N/2 ? j : j-N);

      // Set modes for 0<k<N/2.
      for(k=1;k<N/2;k++)
      {
        pz=k;

        if(cutoff_inigf)
          cutoff = cutoff_process_initial(px,py,pz);
        else
          cutoff = false;

        keff=sqrt(pw2(sin(2.*pi*(double)px/(double)N))+pw2(sin(2.*pi*(double)py/(double)N))+pw2(sin(2.*pi*(double)pz/(double)N)))/dx;
        growAmp >> grow_amp;
        growFac >> grow_fac;
        stabAmp >> stab_amp;
        stabFac >> stab_fac;
        set_mode_gf(px,py,pz,keff,grow_amp,grow_fac,stab_amp,stab_fac,i,j,k,0,cutoff);
      }

      // Set modes with k=0 or N/2.
      // The complex conjugates of these modes appear explicitly on the lattice and must be set to satisfy f(-p)=f(p)*
      if(j>N/2 || (i>N/2 && (j==0 || j==N/2)))
      {
        jconj=(j==0 ? 0 : N-j);

        // k=0
        pz = 0.;
        if(cutoff_inigf)
          cutoff = cutoff_process_initial(px,py,pz);
        else
          cutoff = false;

        keff=sqrt(pw2(sin(2.*pi*(double)px/(double)N))+pw2(sin(2.*pi*(double)py/(double)N)))/dx;
        growAmp >> grow_amp;
        growFac >> grow_fac;
        stabAmp >> stab_amp;
        stabFac >> stab_fac;

        set_mode_gf(px,py,0.,keff,grow_amp,grow_fac,stab_amp,stab_fac,i,j,0,0,cutoff);
        for(fld=2;fld<=4;fld++)
        {
          f[fld][iconj][jconj][0]=f[fld][i][j][0];
          f[fld][iconj][jconj][1]=-f[fld][i][j][1];
          fd[fld][iconj][jconj][0]=fd[fld][i][j][0];
          fd[fld][iconj][jconj][1]=-fd[fld][i][j][1];
        }

        // k=N/2
        pz = (double)N/2.;
        if(cutoff_inigf)
          cutoff = cutoff_process_initial(px,py,pz);
        else
          cutoff = false;
        
        keff=sqrt(pw2(sin(2.*pi*(double)px/(double)N))+pw2(sin(2.*pi*(double)py/(double)N)))/dx;
        growAmp >> grow_amp;
        growFac >> grow_fac;
        stabAmp >> stab_amp;
        stabFac >> stab_fac;

        set_mode_gf_nyquist(px,py,keff,grow_amp,grow_fac,stab_amp,stab_fac,i,j,0,cutoff);
        for(fld=0;fld<=2;fld++)
        {
          fnyquist_gf[fld][iconj][2*jconj]=fnyquist_gf[fld][i][2*j];
          fnyquist_gf[fld][iconj][2*jconj+1]=-fnyquist_gf[fld][i][2*j+1];
          fdnyquist_gf[fld][iconj][2*jconj]=fdnyquist_gf[fld][i][2*j];
          fdnyquist_gf[fld][iconj][2*jconj+1]=-fdnyquist_gf[fld][i][2*j+1];
        }
      }
      else if((i==0 || i==N/2) && (j==0 || j==N/2)) // The 8 "corners" of the lattice are set to real values
      {
        cutoff = false;
        //k=0
        keff=0.;
        growAmp >> grow_amp;
        growFac >> grow_fac;
        stabAmp >> stab_amp;
        stabFac >> stab_fac;
        if(grow_amp!=0. || grow_fac!=0.)
          printf("ERROR: Contents of coulomb files are wrong.\n");
        if(i==0 && j==0)
        {
          set_mode_gf(0.,0.,0.,0.,grow_amp,grow_fac,stab_amp,stab_fac,0,0,0,0,cutoff);
        }
        else
        {
          set_mode_gf(px,py,0.,0.,grow_amp,grow_fac,stab_amp,stab_fac,i,j,0,1,cutoff);
        }

        //k=N/2
        keff=0.;
        growAmp >> grow_amp;
        growFac >> grow_fac;
        stabAmp >> stab_amp;
        stabFac >> stab_fac;
        if(grow_amp!=0. || grow_fac!=0.)
          printf("ERROR: Contents of coulomb files are wrong.\n");
        set_mode_gf_nyquist(px,py,keff,grow_amp,grow_fac,stab_amp,stab_fac,i,j,0,cutoff);
      }
    }
  }

  // DEBUG: whether all cwf values were used
  if(growAmp.eof() && growFac.eof() && stabAmp.eof() && stabFac.eof())
    printf("OK: Coulomb files were exhausted\n");
  else
  {
    printf("ERROR: Coulomb files were not exhausted\n");
    exit(1);
  }


  // Inverse FFT
  if(seed>=0)
  {
    fftrn((double *)f[2],(double *)fnyquist_gf[0],NDIMS,arraysize,-1);
    fftrn((double *)fd[2],(double *)fdnyquist_gf[0],NDIMS,arraysize,-1);
    fftrn((double *)f[3],(double *)fnyquist_gf[1],NDIMS,arraysize,-1);
    fftrn((double *)fd[3],(double *)fdnyquist_gf[1],NDIMS,arraysize,-1);
    fftrn((double *)f[4],(double *)fnyquist_gf[2],NDIMS,arraysize,-1);
    fftrn((double *)fd[4],(double *)fdnyquist_gf[2],NDIMS,arraysize,-1);
  }


  /* A0 is 0 in initialcondition */
  /*---------------------------------------------------------*/
  LOOP
  {
    f[1][i][j][k]  = 0.;
    fd[1][i][j][k] = 0.;
  }

  /* Post processing */
  /*---------------------------------------------------------*/
  save(0);
  output_parameters();
  printf("OK: Finished initial condition of Gauge fields\n");
  return;
}




