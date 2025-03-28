/*
LATTICEEASY consists of the C++ files ``latticeeasy.cpp,''
``initialize.cpp,'' ``evolution.cpp,'' ``output.cpp,''
``latticeeasy.h,'' ``parameters.h,''. (The distribution also includes
the file ffteasy.cpp but this file is distributed separately and
therefore not considered part of the LATTICEEASY distribution in what
follows.) LATTICEEASY is free. We are not in any way, shape, or form
expecting to make money off of these routines. We wrote them for the
sake of doing good science and we're putting them out on the Internet
in case other people might find them useful. Feel free to download
them, incorporate them into your code, modify them, translate the
comment lines into Swahili, or whatever else you want. What we do want
is the following:
1) Leave this notice (i.e. this entire paragraph beginning with
``LATTICEEASY consists of...'' and ending with our email addresses) in
with the code wherever you put it. Even if you're just using it
in-house in your department, business, or wherever else we would like
these credits to remain with it. This is partly so that people can...
2) Give us feedback. Did LATTICEEASY work great for you and help
your work?  Did you hate it? Did you find a way to improve it, or
translate it into another programming language? Whatever the case
might be, we would love to hear about it. Please let us know at the
email address below.
3) Finally, insofar as we have the legal right to do so we forbid
you to make money off of this code without our consent. In other words
if you want to publish these functions in a book or bundle them into
commercial software or anything like that contact us about it
first. We'll probably say yes, but we would like to reserve that
right.

For any comments or questions you can reach us at
gfelder@email.smith.edu
Igor.Tkachev@cern.ch
*/

#include "latticeeasy.hpp"


double f[nflds][N][N][N],fd[nflds][N][N][N];            // Phi,A_0,A_1,A_2,A_3.
double fnyquist_gf[3][N][2*N],fdnyquist_gf[3][N][2*N];  // These are 'gauge mode value strage' for initialize_gf().

double t,t0;

double a=1.,ad=0.;
double asub1=0.,adsub1=0.;
double asub2=0.,adsub2=0.;
double ka=0.,kad=0.;

double hubble_init=0.;
int run_number; 
char mode_[10]="w";
char ext_[500]="_0.dat";
// false->no exit true->exit
bool exitpara = false;


// PBH mass algorithm
/*
double peff_aveList[(int)(1.73205*(N/2))+1];  // thesis method
double pave_effList[(int)(1.73205*(N/2))+1];  // original method
double *pList[2] = {peff_aveList,pave_effList};
double NList[noutput_times+5];
double HList[noutput_times+5];
double aHList[noutput_times+5];*/


int main()
{
  int numsteps=0,output_interval;
  int update_time; 

  /* DEBBUG : seed and thread */
  if(seed<1)
    printf("ERROR: The parameter seed has been set to %d. For correct output set seed to a positive integer.\n",seed);
  if(num_thread>omp_get_max_threads())
  {
    printf("ERROR: Number of threads should below %d.\n",omp_get_max_threads());
    exit(1);
  }

  /* Initialization */
  initialize_ax();
  initialize_gf();
  t=t0;

  output_interval = (int)((tf-t0)/dt)/noutput_times + 1;
  update_time=time(NULL)+print_interval;

  double (*fsub1)[N][N][N]  = (double(*)[N][N][N])malloc(sizeof(double)*nflds*N*N*N);
  double (*fdsub1)[N][N][N] = (double(*)[N][N][N])malloc(sizeof(double)*nflds*N*N*N);
  double (*fsub2)[N][N][N]  = (double(*)[N][N][N])malloc(sizeof(double)*nflds*N*N*N);
  double (*fdsub2)[N][N][N] = (double(*)[N][N][N])malloc(sizeof(double)*nflds*N*N*N);
  double (*kf)[N][N][N]     = (double(*)[N][N][N])malloc(sizeof(double)*nflds*N*N*N);
  double (*kfd)[N][N][N]    = (double(*)[N][N][N])malloc(sizeof(double)*nflds*N*N*N);

  /* Time Evolution */
  while(t<=tf)
  {
    evolve_rk4(dt,SUB,kf,kfd);

    /* Save Values */
    numsteps++;
    if(numsteps%output_interval == 0 && t<tf)
      save(0);

    /* Command Line Output */
    if(time(NULL) >= update_time)
    {
      if(screen_updates)
        printf("%f\n",t);
      update_time += print_interval;
    }

    /* Free memory with error */
    if(exitpara)
    {
      FREE;
      output_parameters();
      printf("ERROR: LATTICEEASY program unfinished\n");
      exit(1);
    }
  }

  /* Post Processing */
  printf("OK: Saving final data\n");
  save(1);
  output_parameters();
  /*
  if(spbh)
  {
    pbh_mass();
  }*/
  FREE;
  printf("OK: LATTICEEASY program finished\n");

  return(0);
}
