#ifndef _PARAMETERS_
#define _PARAMETERS_

// 0->chaotic  1->step
#define MODEL 0

// Model specific adjustable parameters - This section should be replaced by the appropriate one for the model being used.
// (These definitions can be left in model.h but by putting them here you have all the adjustable parameters in one file.)
// ---Adjustable parameters--- //
const double m =0.51e-5;

// Adjustable run parameters
#define NDIMS 3
const int N = 256;
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
const int noutput_times=1000; // Number of times at which to calculate and save output variables
const int print_interval=1; // Interval in seconds between successive outputs of time
const int screen_updates=1; // Set to 1 for time to be periodically output to screen (0 otherwise)
const double checkpoint_interval=0.1; // How often to output a grid image and perform infrequent calculations (see list below). Only done at end if checkpoint_interval=0.
const double store_lattice_times[]={0.}; // An optional list of times at which to close the grid image file and open a new one.
// The variables s<name> control what will be saved (1=save, 0=don't save)
const int smeansvars=1; // Output means and variances. This function is also used to check for exponential instability, so it is generally wise to leave it on.
const int sexpansion=1; // Output scale factor, Hubble constant, and a'' (This is ignored except in self-consistent expansion.)
// The following calculations are performed at intervals given by checkpoint_interval
// 背景量の再現には0で良い
const double t_start_output=0.; // Time to start doing these calculations. This can be reset for any individual calculation.
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