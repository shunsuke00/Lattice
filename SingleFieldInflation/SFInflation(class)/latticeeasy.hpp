/*
This file contains declarations of global variable and important function.
Important function is used in "latticeeasy.cpp".
*/

#ifndef _LATTICEEASYHEADER_
#define _LATTICEEASYHEADER_

#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <float.h>
#include <omp.h>
#include <map>
#include <iostream>

const double pi = (double)(2.*asin(1.));
inline double pw2(double x) {return x*x;}

/////////////////////////////////INCLUDE ADJUSTABLE PARAMETERS///////////////////
#include "parameters.hpp"

/////////////////////////////////GLOBAL DYNAMIC VARIABLES////////////////////////
extern double t,t0;
extern double a,ad;

extern double hubble_init;
extern int run_number;
extern double rescaling;
extern char ext_[500];
extern char mode_[]; // Mode in which to open files, i.e. write ("w") or append ("a+").
extern bool exitpara;

/////////////////////////////////CLASS DEFINITION////////////////////////
class FieldBase{
  protected:
    double (*listall)[nflds][N][N][N];
    double (*f)[N][N][N];
    double (*fd)[N][N][N];

  public:
    FieldBase(double flist[8][nflds][N][N][N])
    {
      listall = flist;
      f=listall[0];
      fd=listall[1];
    }

    // model.hpp
    double potential_energy();
    double dvdf(int fld, int i, int j, int k);

    // Initialize.cpp
    friend void set_mode(double keff, double *field, double *deriv, int real);
    void initialize();

    // element_out.hpp
    double lapl_out(int fld, int i, int j, int k);
    double laplb_out(int fld, int i, int j, int k);
    double field_value(int fld);
    double field_deriv(int fld);
    double field_value2(int fld);
    double field_deriv2(int fld);
    double field_lapl(int fld);
    double fieldderiv_lapl(int fld);
    double fieldderiv_dvdf(int fld);
    // Output.cpp
    void meansvars(int flush);
    void scale(int flush);
    void energy();
    void spectra();
    void pdf();
    void save(int force);
};

class Field:public FieldBase
{
  protected:
    double (*fsub1)[N][N][N];
    double (*fdsub1)[N][N][N];
    double (*fsub2)[N][N][N];
    double (*fdsub2)[N][N][N];
    double (*kf)[N][N][N];
    double (*kfd)[N][N][N];

    double asub1,adsub1;
    double asub2,adsub2;
    double ka,kad;


  public:
    Field(double flist[8][nflds][N][N][N]):FieldBase(flist)
    {
      fsub1 = listall[2];
      fdsub1 = listall[3];
      fsub2 = listall[4];
      fdsub2 = listall[5];
      kf = listall[6];
      kfd = listall[7];

      asub1=0.,adsub1=0.;
      asub2=0.,adsub2=0.;
      ka=0.,kad=0.;
    }

    // model.hpp(overload)
    double potential_energy(int kloop);
    double dvdf(int fld, int i, int j, int k, int kloop);

    // element_evo.hpp
    double lapl(int fld, int i, int j, int k, int kloop_arg);
    double laplb(int fld, int i, int j, int k, int kloop_arg);
    double gradient_energy(int fld, int kloop_arg);
    double stress_3pressure(int kloop_arg);
    // evolution.cpp
    void evolve_fields(double d, int kloop_arg);
    void evolve_derivs(double d, int kloop_arg);
    void evolve_scale(double d, int kloop_arg);
    void evolve_scaled(double d, int kloop_arg);
    void evolve_rk4(double dt);
};

/////////////////////////////////NON-ADJUSTABLE VARIABLES////////////////////////
const double dx=L/(double)N;

/////////////////////////////////DIMENSIONAL SPECIFICATIONS//////////////////////
const int gridsize=N*N*N;

#define FIELD(fld) f[fld][i][j][k]
#define FIELDD(fld) fd[fld][i][j][k]
#define FIELDSUB1(fld) fsub1[fld][i][j][k]
#define FIELDDSUB1(fld) fdsub1[fld][i][j][k]
#define FIELDSUB2(fld) fsub2[fld][i][j][k]
#define FIELDDSUB2(fld) fdsub2[fld][i][j][k]
#define KFIELD(fld) kf[fld][i][j][k]
#define KFIELDD(fld) kfd[fld][i][j][k]

#define LOOP for(i=0;i<N;i++) for(j=0;j<N;j++) for(k=0;k<N;k++)
#define INDEXLIST int i, int j, int k
#define DECLARE_INDICES int i,j,k;
#define FREE do { \
  printf("OK: Freeing memory\n"); \
  free(flist); \
} while(0)


/////////////////////////////////INCLUDE SPECIFIC MODEL//////////////////////////
#include "model.hpp"

/////////////////////////////////FUNCTION DECLARATIONS///////////////////////////
// output.cpp
void output_parameters();


#endif

