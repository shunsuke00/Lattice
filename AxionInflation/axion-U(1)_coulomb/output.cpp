/* This file implements Output function "output_parameters()" and "save()" and "pbh_mass()" */

#include "latticeeasy.hpp"
#include "element_out.hpp"
#include "spectra.hpp"

//char name_[500]; // declare in spectra.hpp

// Output Mean Values Related to Inflaton(Axion) Field
void background(int flush)
{
  static FILE *means_;
  DECLARE_INDICES
  double ave,var,ske,kur;
  double ave_fd;

  static int first=1;
  if(first)
  {
    sprintf(name_,"background%s",ext_);
    means_=fopen(name_,mode_);
    first=0;
  }


  fprintf(means_,"%f",t);

  /* Mean values */
  ave    = field_value(0);
  ave_fd = field_deriv(0);

  /* cumulant (2nd ~ 4th) */
  var=0.,ske=0.,kur=0.;
  #pragma omp parallel for reduction(+:var,ske,kur) num_threads(num_thread) private(j,k)
  LOOP
  {
    var += pw2(AXION-ave);
    ske += pow(AXION-ave,3.);
    kur += pow(AXION-ave,4.);
  }
  var = var/(double)gridsize;
  ske = ske/(double)gridsize;
  kur = kur/(double)gridsize;

  /* Output values */
  fprintf(means_," %e %e %e %e %e\n",ave,ave_fd*pow(a,rescale_s-1.),var,ske/pow(var,3./2.),kur/pw2(var)-3.);
  if(flush)
    fflush(means_);

  /* Check for instability. See whether the field grows exponentially. */
  if(ave+DBL_MAX==ave || (ave!=0. && ave/ave!=1.))
  {
    printf("ERROR: Unstable solution developed. Inflaton field is not numerical at t=%f\n",t);
    fflush(means_);
    exitpara = true;
  }
}


// Outputs Values related to Scale Factor
void scale(int flush)
{
  static FILE *sf_;

  static int first=1;
  if(first)
  {
    sprintf(name_,"sf%s",ext_);
    sf_=fopen(name_,mode_);
    first=0;
  }

  fprintf(sf_,"%f %f %e %f\n",t,a,ad*pow(a,rescale_s-2.),log(a));
  if(flush)
    fflush(sf_);

  // PBH mass algorithm
  /*
  if(spbh)
  {
    static int index = 0;
    NList[index] = log(a);
    HList[index] = ad*pow(a,rescale_s-2.);
    aHList[index] = ad*pow(a,rescale_s-1.);
    index += 1;
  }*/
}


// Output Energy of axion and gauge field + Output Mean Values Related to Backreaction
void energy(int flush)
{
  static FILE *energy_,*coupling_;
  DECLARE_INDICES
  double deriv_energy=0.,grad_energy=0.,pot_energy=0.;
  double axion_energy=0.,axion_pressure=0.,gauge_energy=0.,total=0.;
  double slow_eps0=0.,slow_eps=0.,slow_eta0=0.,slow_eta=0.;

  double ave_fd=0.;
  double xi,cond1,cond2;

  static int first=1;
  if(first)
  {
    sprintf(name_,"energy%s",ext_);
    energy_=fopen(name_,mode_);
    sprintf(name_,"coupling%s",ext_);
    coupling_=fopen(name_,mode_);
    first=0;
  }

  /* energy_0.dat */
  /*-------------------------------------------------------*/
  fprintf(energy_,"%f",t);

  // Energy of Inflaton(Axion)
  deriv_energy   = 0.5*field_deriv2(0)/pw2(a);
  grad_energy    = 0.5*field_lapl(0)/pow(a,2.*rescale_s+2.);
  pot_energy     = potential_energy();
  axion_energy   = deriv_energy + grad_energy + pot_energy;
  axion_pressure = deriv_energy - grad_energy/3. - pot_energy;
  fprintf(energy_," %e %e %e",deriv_energy,grad_energy,pot_energy);

  // Energy of Gauge Field
  gauge_energy = energy_gf();
  fprintf(energy_," %e",gauge_energy);

  // Total and Consistency
  total = axion_energy + gauge_energy;
  fprintf(energy_," %e %e",total,3.*pw2(ad/pw2(a))/total-1.);

  // Slow-roll parameters
  slow_eps0 = (deriv_energy-pot_energy)/(deriv_energy+pot_energy);
  slow_eps0 = 3.*(1.+slow_eps0)/2.;
  slow_eta0 = pow(a,3.)*fieldderiv_dvdf(0)/(ad*deriv_energy*2.*pw2(a));
  slow_eta0 = -6.+2.*slow_eps0-2.*slow_eta0;

  slow_eps = (axion_pressure+gauge_energy/3.)/total;
  slow_eps = 3.*(1.+slow_eps)/2.;
  slow_eta = eta(gauge_energy,slow_eps);
  fprintf(energy_," %e %e %e %e",slow_eps0,slow_eps,slow_eta0,slow_eta);

  // Gauge condition
  if(sgauge)
    fprintf(energy_," %e",gauge_condition());
  
  fprintf(energy_,"\n");
  if(flush)
    fflush(energy_);

  /* coupling_0.dat */
  /*-------------------------------------------------------*/
  fprintf(coupling_,"%f",t);

  // effective coupling
  ave_fd = field_deriv(0);
  xi = g*ave_fd*a/(2.*ad);

  // Backreaction condition (Approximation)
  cond1 = m*pow(a,rescale_s-3.)*pw2(ad)*exp(pi*xi)/(26.*pi*ave_fd*pow(xi,3./2.));
  cond2 = m*pow(a,rescale_s-2.)*ad*exp(pi*xi)/(146.*pow(xi,3./2.));

  fprintf(coupling_," %e %e %e",xi,cond1,cond2);

  // Backreaction condition (Lattice Mean)
  cond1 = backreaction_condition();
  cond2 = pow(a,4.)*gauge_energy/(3.*pw2(ad));
  fprintf(coupling_," %e %e\n",cond1,cond2);

  if(flush)
    fflush(coupling_);
}



// Output Probability Distribution Function of Axion
void pdf_ax()
{
  static FILE *pdf_;

  // bin spacing 0.1sigma
  const int numbins = pdfrange*2*10;
  int numpoints[numbins];
  double valuelist[numbins];

  double sigma,ave;
  double value;

  int i,j,k;

  static int first=1;
  if(first)
  {
    sprintf(name_,"pdf_ax%s",ext_);
    pdf_=fopen(name_,mode_);
    first = 0;
  }

  for(i=0;i<numbins;i++)
  {
    valuelist[i]=i*0.1-0.5*numbins*0.1;
  }


  for(i=0;i<numbins;i++)
  {
    numpoints[i]=0;
  }

  ave = field_value(0);

  sigma=0.;
  LOOP
    sigma += pw2(AXION-ave);
  sigma = sqrt(sigma/(double)gridsize);

  LOOP
  {
    value = (AXION-ave)/sigma;
    if(abs(value)<pdfrange)
    {
      numpoints[(int)(value*10+0.5*numbins)] += 1;
    }
  }

  for(i=0;i<numbins;i++)
  {
    fprintf(pdf_,"%f %d %e\n",valuelist[i],numpoints[i],(double)numpoints[i]/(double)gridsize);
  }
  fprintf(pdf_,"\n");
  fflush(pdf_);
}



/*
// Output an image of all fields and derivatives on the lattice (and a few other variables) to a binary file
void checkpoint()
{
  static FILE *grid_; // File for storing lattice images
  // The following variables are all used for keeping track of when to close one field value file and begin writing to another.
  static int numtimes=sizeof(store_lattice_times)/sizeof(double); // Number of times at which to switch grid image files
  static int current=0; // Index of time value at which to switch files
  static int open; // Indicates whether grid image file is open. (It should be unless the program failed to open the file.)
  int itime; // Integer value of the time up to which a file will be used. Used for generating file names.
  char filename[500];

  static int first=1;
  if(first)
  {
    if(numtimes>0 && store_lattice_times[0]==0.)
      numtimes=0;
    if(numtimes==0) // If no intermediate times are being output simply call the file grid.img
      grid_=fopen("grid.img","wb");
    else // Otherwise label file by the final time it will record data at.
    {
      sprintf(filename,"grid%d.img",(int)store_lattice_times[0]);
      grid_=fopen(filename,"wb");
    }
    first=0;
  }
  else if(current<numtimes && t>store_lattice_times[current]) // If one of the times listed above has been passed switch to another field value file.
  {
    if(open)
      fclose(grid_);
    current++;
    itime = (current<numtimes ? (int)store_lattice_times[current] : (int)tf); // After last time indicated name file with final time of run
    sprintf(filename,"grid%d.img",itime); // Files are labeled by the final time when they will be used, not by the time when they are opened
    grid_=fopen(filename,"wb");
  }

  if(grid_==NULL)
  {
    printf("Error: Grid checkpointing file not open\n");
    open=0;
  }
  else
    open=1;

  if(open) // Write a binary image of all current values to the file grid_.
  {
    rewind(grid_);
    fwrite(&run_number,sizeof(run_number),1,grid_);
    fwrite(&t,sizeof(t),1,grid_);
    fwrite(&a,sizeof(a),1,grid_);
    fwrite(&ad,sizeof(ad),1,grid_);
    fwrite(f,sizeof(double),nflds*gridsize,grid_);
    fwrite(fd,sizeof(double),nflds*gridsize,grid_);
    fflush(grid_);
  }
}
*/

// Convert a time in seconds to a more readable form and print the results
void readable_time(int t, FILE *info_)
{
  int tminutes=60,thours=60*tminutes,tdays=24*thours;

  if(t==0)
  {
    fprintf(info_,"less than 1 second\n");
    return;
  }

  // Days
  if(t>tdays)
  {
    fprintf(info_,"%d days",t/tdays);
    t = t%tdays;
    if(t>0)
      fprintf(info_,", ");
  }  
  // Hours
  if(t>thours)
  {
    fprintf(info_,"%d hours",t/thours);
    t = t%thours;
    if(t>0)
      fprintf(info_,", ");
  }
  // Minutes
  if(t>tminutes)
  {
    fprintf(info_,"%d minutes",t/tminutes);
    t = t%tminutes;
    if(t>0)
      fprintf(info_,", ");
  }
  // Seconds
  if(t>0)
    fprintf(info_,"%d seconds",t);
  fprintf(info_,"\n");
  return;
}

/////////////////////////////////////////////////////
// Important Function
/////////////////////////////////////////////////////

// Output information about the run parameters.
void output_parameters()
{
  static FILE *info_,*python_;
  static time_t tStart,tFinish;

  static int first=1;
  if(first) // At beginning of run
  {
    /* info_0.dat */
    sprintf(name_,"info%s",ext_);
    info_=fopen(name_,mode_);

    fprintf(info_,"-----------------------------\n");
    fprintf(info_,"Axion - U(1) Model\n");
    fprintf(info_,"-----------------------------\n");
    fprintf(info_,"coupling alpha/f = %f\n",g);

    fprintf(info_,"\n-----------------------------\n");
    fprintf(info_,"Inflaton Potential\n");
    fprintf(info_,"-----------------------------\n");
    modelinfo(info_);

    fprintf(info_,"\n-----------------------------\n");
    fprintf(info_,"General Program Information\n");
    fprintf(info_,"-----------------------------\n");
    fprintf(info_,"Grid size=%d^%d\n",N,NDIMS);
    fprintf(info_,"L=%f\n",L);
    fprintf(info_,"dt=%f, dt/dx=%f\n",dt,dt/dx);
    fprintf(info_,"\n-----------------------------\n");
    fprintf(info_,"Detailed Program Setup\n");
    fprintf(info_,"-----------------------------\n");
    if(hankel)
      fprintf(info_,"Hankel Fluctuation Initial Condition (Axion)\n");
    else
      fprintf(info_,"EXP Fluctuation Initial Condition (Axion)\n");
    if(coulomb_po)
    {
      if(coulomb_ne)
        fprintf(info_,"Coulomb Wave Function (Both Polarization) Initial Condition (Gauge)\n");
      else
        fprintf(info_,"Coulomb Wave Function (Positive Polarization) Initial Condition (Gauge)\n");
    }
    else
    {
      if(coulomb_ne)
        fprintf(info_,"Coulomb Wave Function (Negative Polarization) Initial Condition (Gauge)\n");
      else
        fprintf(info_,"EXP Initial Condition (Gauge)\n");
    }
    fprintf(info_,"Polarization=%d\n",POLARIZATION);
    fprintf(info_,"Cutoff method=%d\n",cutoff_method);
    fprintf(info_,"initial cutoff=(%d,%d)\n",cutoff_iniax,cutoff_inigf);
    fprintf(info_,"spectra cutoff=(%d,%d)\n",cutoff_speax,cutoff_spegf);
    fprintf(info_,"seed=%d\n",seed);
    fprintf(info_,"sbackground=%d\n",sbackground);
    fprintf(info_,"sexpansion=%d\n",sexpansion);
    fprintf(info_,"senergy=%d\n",senergy);
    fprintf(info_,"sgauge=%d\n",sgauge);
    fprintf(info_,"sspectra=%d\n",sspectra);
    fprintf(info_,"spdf=%d\n",spdf);
    if(spdf)
      fprintf(info_,"pdfrange=%d\n",pdfrange);
    fprintf(info_,"\n-----------------------------\n");
    /*fprintf(info_,"Input parameter for PBH mass\n");
    fprintf(info_,"-----------------------------\n");
    fprintf(info_,"spbh=%d\n",spbh);
    if(spbh)
    {
      fprintf(info_,"Hcmb/m=%f\n",Hcmb);
      fprintf(info_,"N0-Ncmb=%f\n",Ncmb);
      fprintf(info_,"kcmb=%f Mpc^-1\n",kcmb);
      fprintf(info_,"effeciency factor gamma=%f\n",gamm);
      fprintf(info_,"effective degree of freedom g_*=%f\n\n",gsta);
    }*/
    time(&tStart);
    fprintf(info_,"\nRun began at %s",ctime(&tStart));

    /* python.dat */
    python_=fopen("python.dat","w");
    fprintf(python_,"%d %d %d %d %f %e %d %d %d %d %d %d %d\n",
    BACKREACTION,POLARIZATION,MODEL,N,L,m,
    sbackground,sexpansion,senergy,sgauge,sspectra,spdf,pdfrange);
    fprintf(python_,"%f %f %f ",g,initfield[0],ad);

    first=0;
  }
  else // At end of run
  {
    /* info_0.dat */
    time(&tFinish);
    fprintf(info_,"Run ended at %s",ctime(&tFinish));
    fprintf(info_,"\nRun from t=%f to t=%f took ",t0,t);
    readable_time((int)(tFinish-tStart),info_);
    fprintf(info_,"\n");

    /* python.dat */
    fprintf(python_,"%f %f %f\n",a,ad,field_deriv(0));
  }

  fflush(info_);
  fflush(python_);
  return;
}

// Calculate and save quantities (means, variances, etc.). If force>0 all infrequent calculations will be performed
void save(int force)
{
  static double tsave=t0; // Used to indicate when output files should be flushed.
  int flush=0; // Indicates whether to flush files, checkpoint, and perform certain infrequent calculations

  if(checkpoint_interval==0.) // In this case only checkpoint at the end of the run
    tsave=tf;

  if(t>=tsave || force>0) // Check if the program should flush
  {
    tsave = t + checkpoint_interval;
    flush=1;
  }

  // Frequent calculation
  if(sbackground)
    background(flush);
  if(sexpansion)
    scale(flush);
  if(senergy)
    energy(flush);

  // Infrequent calculations
  if(sspectra && flush && t>=tspectra)
  {
    spectra_ax();
    spectra_gf();
    if(spdf)
      pdf_ax();
  }
  
  /*
  if(scheckpoint && flush && t>=tcheckpoint)
    checkpoint();
  */
}


/*
// PBH mass algorithm
void pbh_mass()
{
  FILE *pbhfiles[2];

  sprintf(name_,"PBHmass_the%s",ext_);
  pbhfiles[0]=fopen(name_,mode_);
  sprintf(name_,"PBHmass_ori%s",ext_);
  pbhfiles[1]=fopen(name_,mode_);

  int index,i;
  double N_be,N_af,H_be,H_af,aH_be,aH_af;
  double deltaN,deltaH;
  double mass,kratio;
  double k;

  for(int n=0;n<2;n++)
  {
    if(n==1)
      i = 0;
    else
      i = 1;
    while(pList[n][i]>0.00001)
    {
      index = 0;
      k = pList[n][i];
      if(k<aHList[0])
      {
        fprintf(pbhfiles[n],"0 %e %f %f %f %f %f %f\n",k,0.,0.,0.,0.,0.,0.);
      }
      else
      {
        while(k>aHList[index])
          index += 1; // indexは境界の後を指す
        N_be = NList[index-1];
        N_af = NList[index];
        H_be = HList[index-1];
        H_af = HList[index];
        aH_be = aHList[index-1];
        aH_af = aHList[index];

        deltaN = (k-aH_be)*(N_af-N_be)/(aH_af-aH_be);
        if(deltaN<=0.)
          fprintf(pbhfiles[n],"1 %e %f %f %f %f %f %f\n",k,0.,0.,0.,0.,0.,0.);
        else
        {
          deltaH = deltaN*(H_af-H_be)/(N_af-N_be);
          if(deltaH>=0.)
            fprintf(pbhfiles[n],"2 %e %f %f %f %f %f %f\n",k,0.,0.,0.,0.,0.,0.);
          else
          {
            mass = 2.79*pow(10.,-12.)*(gamm/0.2)*pow(106.75/gsta,1./6.)*pw2(pow(10.,12.)*Hcmb/(kcmb*(H_be+deltaH)))*exp(-2.*(Ncmb+N_be+deltaN));
            kratio = exp(Ncmb+N_be+deltaN)*(H_be+deltaH)/Hcmb;
            fprintf(pbhfiles[n],"3 %e %f %f %e %e %e %e\n",k,N_be+deltaN,H_be+deltaH,mass,kratio,(N_af-N_be)/N_af,(H_af-H_be)/H_af);
          }
        }
      }
      fflush(pbhfiles[n]);
      i += 1;
    }
  }
}*/