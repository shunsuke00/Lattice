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


double t,t0;
double a=1.,ad=0.;

double hubble_init=0.;
int run_number;
char mode_[10]="w";
double rescaling=1.;
char ext_[500]="_0.dat";
// false->no exit true->exit
bool exitpara = false;


int main()
{
  int numsteps=0,output_interval;
  int update_time;

  /* DEBUG: seed */
  if(seed<1)
    printf("ERROR: The parameter seed has been set to %d. For correct output set seed to a positive integer.",seed);

  /* Dynamical Allocation of Field Values Memory */
  double (*flist)[nflds][N][N][N] = (double(*)[nflds][N][N][N])malloc(8*sizeof(double)*nflds*N*N*N);
  Field field(flist);

  /* Initialization */
  field.initialize();
  t=t0;

  output_interval = (int)((tf-t0)/dt)/noutput_times + 1;
  update_time=time(NULL)+print_interval;

  /* Evolution */
  while(t<=tf)
  {
    field.evolve_rk4(dt);

    /* Save values */
    numsteps++;
    if(numsteps%output_interval == 0 && t<tf)
      field.save(0);

    /* Command line output */
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
  field.save(1);
  output_parameters();
  FREE;
  printf("OK: LATTICEEASY program finished\n");

  return(0);
}
