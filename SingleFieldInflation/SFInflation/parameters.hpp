#ifndef _PARAMETERS_
#define _PARAMETERS_


#define MODEL 0           // Inflaton potential choice (0->chaotic  1->step)
const double m =0.51e-5;  // Inflaton mass[M_pl]

// Adjustable run parameters
#define NDIMS 3
const int N = 64;
const int nflds = 1;
#if MODEL==0
const double L = 1.4;
const double initfield[]={14.5};
const double initderivs[]={-0.8152128};
#elif MODEL==1
const double L = 0.6;
const double initfield[]={14.5};
const double initderivs[]={-0.8233865};
#endif
const double dt = .0001;
const double tf=1.2;
const int seed=2;

const bool hankel=true; // Inflaton fluctuation initial setup
const double kcutoff=0; // Momentum for initial lowpass filter. Set to 0 to not filter

// If and how to continue previous runs.
// If no grid image is available in the run directory then a new run will be started irrespective of continue_run.
// 0=Start a new run at t=0. 1=Continue old run, appending new data to old data files. 2=Continue old run, creating new data files for all output. (Old ones will not in general be overwritten.)
const int continue_run=0;

// Variables controlling output
const char alt_extension[]="";
const int noutput_times=1000;
const int print_interval=1;
const int screen_updates=1;
const double checkpoint_interval=0.1;
const double store_lattice_times[]={0.}; // An optional list of times at which to close the grid image file and open a new one.
// The variables s<name> control what will be saved (1=save, 0=don't save)
const int smeansvars=1;
const int sexpansion=1;
// The following calculations are performed at intervals given by checkpoint_interval
const double t_start_output=0.;
const int senergy=1;
  const double tenergy=t_start_output;
const int sspectra=1;
  const double tspectra=t_start_output;
  const int spdf=1;
  const int pdfrange=5;

/*
const int scheckpoint=0;
  const double tcheckpoint=t_start_output;
*/

// The number of threads in openMP
const int num_thread=8;

#endif