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
#include <fstream>
#include <iomanip>

const double pi = (double)(2.*asin(1.));
inline double pw2(double x) {return x*x;}

////////////////INCLUDE ADJUSTABLE PARAMETERS AND POLARIZATION VECTOR/////////////
#include "parameters.hpp"
#include "polarization.hpp"

/////////////////////////////////GLOBAL DYNAMIC VARIABLES////////////////////////
extern double t,t0;

extern double a,ad;
extern double asub1,adsub1;
extern double asub2,adsub2;
extern double ka,kad;

extern double hubble_init;
extern int run_number;
extern char ext_[500];
extern char mode_[]; // Mode in which to open files, i.e. write ("w") or append ("a+").
extern bool exitpara;


// PBH mass algorithm
/*
extern double peff_aveList[(int)(1.73205*(N/2))+1];  // thesis method
extern double pave_effList[(int)(1.73205*(N/2))+1];  // original method
extern double *pList[2];
extern double NList[noutput_times+5];
extern double HList[noutput_times+5];
extern double aHList[noutput_times+5];*/


/////////////////////////////////NON-ADJUSTABLE VARIABLES////////////////////////
//const double dx=L/(double)N;  //declare in polarization.hpp

/////////////////////////////////DIMENSIONAL SPECIFICATIONS//////////////////////
extern double f[nflds][N][N][N],fd[nflds][N][N][N];
extern double fnyquist_gf[3][N][2*N],fdnyquist_gf[3][N][2*N];

const int gridsize=N*N*N;

/////////////////////////////////MACRO///////////////////////////////////////////
/* Field macro */
// fld=0->Inflaton(axion)
#define AXION      f[0][i][j][k]
#define AXIOND     fd[0][i][j][k]
#define AXIONSUB1  fsub1[0][i][j][k]
#define AXIONDSUB1 fdsub1[0][i][j][k]
#define AXIONSUB2  fsub2[0][i][j][k]
#define AXIONDSUB2 fdsub2[0][i][j][k]
#define KAXION     kf[0][i][j][k]
#define KAXIOND    kfd[0][i][j][k]
// fld=1~4->Gauge, com=0~3
#define GAUGE(com)      f[com+1][i][j][k]
#define GAUGED(com)     fd[com+1][i][j][k]
#define GAUGESUB1(com)  fsub1[com+1][i][j][k]
#define GAUGEDSUB1(com) fdsub1[com+1][i][j][k]
#define GAUGESUB2(com)  fsub2[com+1][i][j][k]
#define GAUGEDSUB2(com) fdsub2[com+1][i][j][k]
#define KGAUGE(com)     kf[com+1][i][j][k]
#define KGAUGED(com)    kfd[com+1][i][j][k]
// whole field
#define FIELD(fld)      f[fld][i][j][k]
#define FIELDD(fld)     fd[fld][i][j][k]
#define FIELDSUB1(fld)  fsub1[fld][i][j][k]
#define FIELDDSUB1(fld) fdsub1[fld][i][j][k]
#define FIELDSUB2(fld)  fsub2[fld][i][j][k]
#define FIELDDSUB2(fld) fdsub2[fld][i][j][k]
#define KFIELD(fld)     kf[fld][i][j][k]
#define KFIELDD(fld)    kfd[fld][i][j][k]

// Argument of functions
#define KSUB_ARG double kf[nflds][N][N][N], double kfd[nflds][N][N][N]
#define SUB_ARG double fsub1[nflds][N][N][N], double fdsub1[nflds][N][N][N], double fsub2[nflds][N][N][N], double fdsub2[nflds][N][N][N]
#define SUB fsub1,fdsub1,fsub2,fdsub2
#define FSUB_ARG double fsub1[nflds][N][N][N], double fsub2[nflds][N][N][N]
#define FSUB fsub1,fsub2

/* Other macro */
#define FIELDPOINT(fld,i,j,k) f[fld][i][j][k]
#define LOOP for(i=0;i<N;i++) for(j=0;j<N;j++) for(k=0;k<N;k++)
#define INDEXLIST int i, int j, int k
#define DECLARE_INDICES int i,j,k;
#define FREE do { \
  printf("OK: Freeing memory\n"); \
  free(fsub1); \
  free(fdsub1); \
  free(fsub2); \
  free(fdsub2); \
  free(kf); \
  free(kfd); \
} while(0)

/*
#define CUTOFF_PROCESS_INITIAL do{ \
  if(cutoff_method==0){ \
    cutoff = (sqrt(pw2(px) + pw2(py) + pw2(pz)) > modecutoff); \
  } else if(cutoff_method==1){ \
    cutoff = (abs(px)>(double)N/4. || abs(py)>(double)N/4. || abs(pz)>(double)N/4.); \
  } else{ \
    printf("ERROR: Value of mutoff_method is wrong.\n"); \
  } \
} while(0)*/



/////////////////////////////////INCLUDE SPECIFIC MODEL//////////////////////////
#include "model.hpp"

/////////////////////////////////FUNCTION DECLARATIONS///////////////////////////
// initialize.cpp
void initialize_ax();
void initialize_gf();
// evolution.cpp
void evolve_rk4(double dt,SUB_ARG,KSUB_ARG);
// output.cpp
void output_parameters();
void save(int force);
//void pbh_mass();
// ffteasy.cpp
void fftr1(double f_arg[], int N_arg, int forward);
void fftrn(double f_arg[], double fnyquist[], int ndims, int size[], int forward);

#endif

