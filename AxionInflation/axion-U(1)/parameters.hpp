#ifndef _PARAMETERS_
#define _PARAMETERS_


#define BACKREACTION 2    // Backreaction Choice (0->Weak 1->Strong 2->BumpyAxionInflaiton)
#define POLARIZATION 1    // Polarization Vector Choice (0->Reciprocal 1->EffectiveMomentum)
#define MODEL 1           // Inflaton Potential Choice (0->Chaotic 1->BumpyAxion)
const double m =0.51e-5;  // Inflaton mass[M_pl]

// Adjustable run parameters
#define NDIMS 3
const int N =256;
const int nflds = 5;

#if BACKREACTION==0
const double L = 2.0;
const double initfield[]={-14.5};
const double initderiv[]={0.8152128};
const double g=42.;

const double dt = .0004;
const double tf=1.2;
#elif BACKREACTION==1
const double L = 0.2;
const double initfield[]={-5.5};
const double initderiv[]={0.80815};
const double g=25.;

const double dt = .0004;
const double tf=7.0;
#elif BACKREACTION==2
const double L = 3.5;
//const double initfield[]={-4.6205366};
//const double initderiv[]={0.092055699};
//const double initfield[]={-4.5929167};
//const double initderiv[]={0.12242202};
//const double initfield[]={-4.5042820};//H=1.93
//const double initderiv[]={0.24034564};
const double initfield[]={-4.4273861};//H=1.92
const double initderiv[]={0.36180591};
//const double initfield[]={-4.3387189};
//const double initderiv[]={0.51552635};
const double g=3.3*6.7;

const double dt = .0004;
const double tf=6.0;
#endif

const int seed=4;
const bool hankel = true;     // Inflaton fluctuation initial setup
const bool coulomb_po = true; // Gauge fluctuation initial setup (positive polarization)
const bool coulomb_ne = true; // Gauge fluctuation initial setup (negative polarization)
const std::string path_base = "/Users/tsuchitashunsuke/Dropbox/lattice/coulomb/CWFvalue/";


// 0->reprocical vector amplitude 1->one dimensional cutoff
const int cutoff_method = 0;
// true->cutoff false->no cutoff
const bool cutoff_iniax = true;
const bool cutoff_inigf = true;
const bool cutoff_speax = true;
const bool cutoff_spegf = true;
const double modecutoff = (double)N/4.;

// If and how to continue previous runs.
// If no grid image is available in the run directory then a new run will be started irrespective of continue_run.
// 0=Start a new run at t=0. 1=Continue old run, appending new data to old data files. 2=Continue old run, creating new data files for all output. (Old ones will not in general be overwritten.)
const int continue_run=0;

// Variables controlling output
const char alt_extension[]="";
const int noutput_times=1000;
const int print_interval=1;
const int screen_updates=1;
#if BACKREACTION==0
const double checkpoint_interval=0.1;
#elif BACKREACTION==1
const double checkpoint_interval=0.3;
#elif BACKREACTION==2
const double checkpoint_interval=0.3;
#endif
const double store_lattice_times[]={0.}; // An optional list of times at which to close the grid image file and open a new one.
// The variables s<name> control what will be saved (1=save, 0=don't save)
const int sbackground=1;
const int sexpansion=1;
const int senergy=1;
  const int sgauge=1;
// The following calculations are performed at intervals given by checkpoint_interval
const double t_start_output=0.;
const int sspectra=1;
  const double tspectra=t_start_output;
  const int spdf=1;
  const int pdfrange=5;
// PBH calculation
//const int spbh=1;

/*
const int scheckpoint=0;
  const double tcheckpoint=t_start_output;
*/

// The number of threads in openMP
const int num_thread=18;

#endif