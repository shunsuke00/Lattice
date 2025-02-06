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
extern double asub1,adsub1;
extern double asub2,adsub2;
extern double ka,kad;

extern double hubble_init;
extern int run_number;
extern double rescaling;
extern char ext_[500];
extern char mode_[]; // Mode in which to open files, i.e. write ("w") or append ("a+").
extern bool exitpara;

/////////////////////////////////NON-ADJUSTABLE VARIABLES////////////////////////
const double dx=L/(double)N;

/////////////////////////////////DIMENSIONAL SPECIFICATIONS//////////////////////
extern double f[nflds][N][N][N],fd[nflds][N][N][N];

const int gridsize=N*N*N;

#define FIELD(fld) f[fld][i][j][k]
#define FIELDD(fld) fd[fld][i][j][k]
#define FIELDSUB1(fld) fsub1[fld][i][j][k]
#define FIELDDSUB1(fld) fdsub1[fld][i][j][k]
#define FIELDSUB2(fld) fsub2[fld][i][j][k]
#define FIELDDSUB2(fld) fdsub2[fld][i][j][k]
#define KFIELD(fld) kf[fld][i][j][k]
#define KFIELDD(fld) kfd[fld][i][j][k]

#define SUB_ARG double fsub1[nflds][N][N][N], double fdsub1[nflds][N][N][N], double fsub2[nflds][N][N][N], double fdsub2[nflds][N][N][N]
#define SUB fsub1,fdsub1,fsub2,fdsub2
#define FSUB_ARG double fsub1[nflds][N][N][N], double fsub2[nflds][N][N][N]
#define FSUB fsub1,fsub2


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


/////////////////////////////////INCLUDE SPECIFIC MODEL//////////////////////////
#include "model.hpp"

/////////////////////////////////FUNCTION DECLARATIONS///////////////////////////
// initialize.cpp
void initialize();
// evolution.cpp
void evolve_rk4(double dt, SUB_ARG, double kf[nflds][N][N][N], double kfd[nflds][N][N][N]);
// output.cpp
void output_parameters();
void save(int force);
// ffteasy.cpp
void fftr1(double f_arg[], int N_arg, int forward);
void fftrn(double f_arg[], double fnyquist[], int ndims, int size[], int forward);


#endif

