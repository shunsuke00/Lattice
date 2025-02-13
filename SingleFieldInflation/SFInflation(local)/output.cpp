/* This file implement Output Function "output_parameters()" and "save()" */

#include "latticeeasy.hpp"
#include "element_out.hpp"

char name_[500]; // Filenames - set differently by each function to open output files

// Output Mean Values Related to Inflaton Field
void meansvars(int flush, ARG)
{
  static FILE *means_[nflds];
  DECLARE_INDICES
  int fld;
  double ave,var,ske,kur;
  double ave_fd;

  double (*f)[N][N][N] = flist[0];

  static int first=1;
  if(first)
  {
    for(fld=0;fld<nflds;fld++)
    {
      sprintf(name_,"means%d%s",fld,ext_);
      means_[fld]=fopen(name_,mode_);
    }
    first=0;
  }

  
  for(fld=0;fld<nflds;fld++)
  {
    fprintf(means_[fld],"%f",t);

    /* Mean values */
    ave    = field_value(fld,flist);
    ave_fd = field_deriv(fld,flist);

    /* cumulant (2nd ~ 4th) */
    var=0.,ske=0.,kur=0.;
    #pragma omp parallel for reduction(-:var,ske,kur) num_threads(num_thread) private(j,k)
    LOOP
    {
      var += pw2(FIELD(fld)-ave);
      ske += pow(FIELD(fld)-ave,3.);
      kur += pow(FIELD(fld)-ave,4.);
    }
    var = var/(double)gridsize;
    ske = ske/(double)gridsize;
    kur = kur/(double)gridsize;

    /* Output values */
    fprintf(means_[fld]," %e %e %e %e %e\n",ave,ave_fd*pow(a,rescale_s-1.),var,ske/pow(var,3./2.),kur/pw2(var)-3.);
    if(flush)
      fflush(means_[fld]);

    /* Check for instability. See whether the field has grown exponentially */
    if(ave+DBL_MAX==ave || (ave!=0. && ave/ave!=1.))
    {
      printf("ERROR: Unstable solution developed. Field %d not numerical at t=%f\n",fld,t);
      exitpara = true;
    }
  }
}


// Output Values related to Scale Factor
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
}


// Output Values related to Energy of Inflaton
void energy(ARG)
{
  static FILE *energy_;
  DECLARE_INDICES
  int fld;
  double deriv_energy,grad_energy,pot_energy,total=0.;
  double fd2,fddvdf,grad,dgrad;
  double slow_eps0=0.,slow_eps=0.,slow_eta0=0.,slow_eta=0.;

  static int first=1;
  if(first)
  {
    sprintf(name_,"energy%s",ext_);
    energy_=fopen(name_,mode_);
    first=0;
  }

  fprintf(energy_,"%f",t);

  // Kinetic energy
  for(fld=0;fld<nflds;fld++)
  {
    deriv_energy = 0.5*field_deriv2(fld,flist)/pw2(a);
    total += deriv_energy;
    fprintf(energy_," %e",deriv_energy);
  }

  /* Gradient energy */
  for(fld=0;fld<nflds;fld++)
  {
    grad_energy = 0.5*field_lapl(fld,flist)/pow(a,2.*rescale_s+2.);
    total += grad_energy;
    fprintf(energy_," %e",grad_energy);
  }

  // Potential energy and total energy
  pot_energy = potential_energy(1,flist);
  total += pot_energy;
  fprintf(energy_," %e %e",pot_energy,total);

  // Energy conservation and slow-roll parameter
  slow_eps0 = (deriv_energy-pot_energy)/(deriv_energy+pot_energy);
  slow_eps0 = 3.*(1.+slow_eps0)/2.;
  fd2       = deriv_energy*2.*pw2(a);
  fddvdf    = fieldderiv_dvdf(0,flist);
  slow_eta0 = -6.+2.*slow_eps0-2.*pow(a,3.)*fddvdf/(ad*fd2);

  slow_eps = (deriv_energy-grad_energy/3.-pot_energy)/total;
  slow_eps = 3.*(1.+slow_eps)/2.;
  grad     = grad_energy*pow(a,2.*rescale_s+2.);
  dgrad    = fieldderiv_lapl(0,flist);
  slow_eta = 2.*slow_eps-3.*(1.+fd2/(fd2-pow(a,-2.*rescale_s)*grad/3.))-2.*pow(a,3.)*(fddvdf)/(ad*(fd2-pow(a,-2.*rescale_s)*grad/3.));
  slow_eta = slow_eta - (grad-4.*a*dgrad/ad)/(3.*pw2(a)*(pow(a,2.*rescale_s-2.)*fd2-grad/(3.*pw2(a))));

  fprintf(energy_," %8.7e %e %e %e %e\n",3.*pw2(ad/pw2(a))/total-1.,slow_eps0,slow_eps,slow_eta0,slow_eta);

  fflush(energy_);
}




// Output Inflaton Power Spectra
void spectra(ARG)
{
  static FILE *spectra_[nflds],*spectra_eff_[nflds],*spectratimes_;
  int fld;
  const int maxnumbins=(int)(1.73205*(N/2))+1; // Number of bins
  int numpoints[maxnumbins],numpoints_eff[maxnumbins];
  double p[maxnumbins],pave[maxnumbins],peff_ave[maxnumbins],f2[maxnumbins];
  double p_eff[maxnumbins],pave_eff[maxnumbins],f2_eff[maxnumbins];

  int numbins=(int)(sqrt((double)NDIMS)*(N/2))+1;
  double pmagnitude; // |m|=Sqrt(mx^2+my^2+mz^2)
  double peffmagnitude;
  double dp=2.*pi/L;
  double fp2;

  // parameter with dimension
  int i,j,k,mx,my,mz;
  double fnyquist[N][2*N];
  int arraysize[]={N,N,N};

  double (*f)[N][N][N] = flist[0];


  static int first=1;
  if(first)
  {
    for(fld=0;fld<nflds;fld++)
    {
      sprintf(name_,"spectra%d%s",fld,ext_);
      spectra_[fld]=fopen(name_,mode_);
      sprintf(name_,"spectra_eff%d%s",fld,ext_);
      spectra_eff_[fld]=fopen(name_,mode_);
    }
    sprintf(name_,"spectratimes%s",ext_);
    spectratimes_=fopen(name_,mode_);
    first=0;
  }

  // Pre Processing
  for(i=0;i<numbins;i++)
  {
    p[i]=dp*i;
    p_eff[i]=dp*i;
  }
    
  for(fld=0;fld<nflds;fld++)
  {
    for(i=0;i<numbins;i++)
    {
      numpoints[i]=0;
      pave[i]=0.;
      peff_ave[i]=0.;
      f2[i]=0.;

      numpoints_eff[i]=0;
      pave_eff[i]=0.;
      f2_eff[i]=0.;
    }

    // FFT
    fftrn((double *)f[fld],(double *)fnyquist,NDIMS,arraysize,1);


    // Calculating Spectra
    for(i=0;i<N;i++)
    {
      mx=(i<=N/2 ? i : i-N);
      for(j=0;j<N;j++)
      {
        my=(j<=N/2 ? j : j-N);

        /* Mode with 0<k<N/2 */
        for(k=1;k<N/2;k++) 
        {
          mz=k;

          pmagnitude=sqrt(pw2(mx)+pw2(my)+pw2(mz));
          peffmagnitude=(double)N*sqrt(pw2(sin(pi*(double)mx/(double)N))+pw2(sin(pi*(double)my/(double)N))+pw2(sin(pi*(double)mz/(double)N)))/pi;
          fp2=pw2(f[fld][i][j][2*k])+pw2(f[fld][i][j][2*k+1]);

          numpoints[(int)pmagnitude] += 2;
          pave[(int)pmagnitude] += 2.*pmagnitude;
          peff_ave[(int)pmagnitude] += 2.*peffmagnitude;
          f2[(int)pmagnitude] += 2.*fp2;
    
          numpoints_eff[(int)peffmagnitude] += 2;
          pave_eff[(int)peffmagnitude] += 2.*peffmagnitude;
          f2_eff[(int)peffmagnitude] += 2.*fp2;
        }

        // Mode with k=0,N/2
        for(k=0;k<=N/2;k+=N/2)
        {
          mz=k;

          pmagnitude=sqrt(pw2(mx)+pw2(my)+pw2(mz));
          peffmagnitude=(double)N*sqrt(pw2(sin(pi*(double)mx/(double)N))+pw2(sin(pi*(double)my/(double)N))+pw2(sin(pi*(double)mz/(double)N)))/pi;
          if(k==0) // Amplitudes of k=0 modes
          {   
            fp2=pw2(f[fld][i][j][0])+pw2(f[fld][i][j][1]);
          }
          else // Modes with k=N/2 are stored in the array fnyquist
          {
            fp2=pw2(fnyquist[i][2*j])+pw2(fnyquist[i][2*j+1]);
          }

          numpoints[(int)pmagnitude]++;
          pave[(int)pmagnitude] += pmagnitude;
          peff_ave[(int)pmagnitude] += peffmagnitude;
          f2[(int)pmagnitude] += fp2;

          if(k==0 && j==0 && i==0)
          {
            continue;
          }
          numpoints_eff[(int)peffmagnitude]++;
          pave_eff[(int)peffmagnitude] += peffmagnitude;
          f2_eff[(int)peffmagnitude] += fp2;
        }
      }
    }
    
    /* spectra_0.dat and spectra_eff_0.dat */
    for(i=0;i<numbins;i++)
    {
      // Convert sums to averages. (numpoints[i] should always be greater than zero.)
      if(numpoints[i]>0)
      {
        pave[i] = pave[i]*dp/numpoints[i];
        peff_ave[i] = peff_ave[i]*dp/numpoints[i];
        f2[i] = f2[i]/numpoints[i];
      }
      if(numpoints_eff[i]>0)
      {
        pave_eff[i] = pave_eff[i]*dp/numpoints_eff[i];
        f2_eff[i] = f2_eff[i]/numpoints_eff[i];
      }
      // Output the momentum, number of points, omega, and calculated spectra for each bin
      fprintf(spectra_[fld],"%e %d %e %e %e\n",
        p[i],numpoints[i],pave[i],peff_ave[i],f2[i]);
      fprintf(spectra_eff_[fld],"%e %d %e %e\n",
        p_eff[i],numpoints_eff[i],pave_eff[i],f2_eff[i]);
    }
    fprintf(spectra_[fld],"\n");
    fflush(spectra_[fld]);
    fprintf(spectra_eff_[fld],"\n");
    fflush(spectra_eff_[fld]);

    // Inverse FFT
    fftrn((double *)f[fld],(double *)fnyquist,NDIMS,arraysize,-1);
  }

  /* spectratimes_0.dat */
  double av_fd=field_deriv(0,flist);

  fprintf(spectratimes_,"%f %f %f\n",t,log(a),ad/(a*av_fd));
  fflush(spectratimes_);

  return;
}


// Output Probability Distribution Function
void pdf(ARG)
{
  static FILE *pdf_[nflds];

  // bin spacing 0.1sigma
  const int numbins = pdfrange*2*10;
  int numpoints[numbins];
  double valuelist[numbins];

  double sigma,ave;
  double value;

  int i,j,k,fld;

  double (*f)[N][N][N] = flist[0];

  static int first=1;
  if(first)
  {
    for(fld=0;fld<nflds;fld++)
    {
      sprintf(name_,"pdf%d%s",fld,ext_);
      pdf_[fld]=fopen(name_,mode_);
    }
    first = 0;
  }

  for(i=0;i<numbins;i++)
  {
    valuelist[i]=i*0.1-0.5*numbins*0.1;
  }


  for(fld=0;fld<nflds;fld++)
  {
    for(i=0;i<numbins;i++)
    {
      numpoints[i]=0;
    }

    ave = field_value(fld,flist);

    sigma=0.;
    LOOP
      sigma += pw2(FIELD(fld)-ave);
    sigma = sqrt(sigma/(double)gridsize);

    LOOP
    {
      value = (FIELD(fld)-ave)/sigma;
      if(abs(value)<pdfrange)
      {
        numpoints[(int)(value*10+0.5*numbins)] += 1;
      }
    }

    for(i=0;i<numbins;i++)
    {
      fprintf(pdf_[fld],"%f %d %e\n",valuelist[i],numpoints[i],(double)numpoints[i]/(double)gridsize);
    }
    fprintf(pdf_[fld],"\n");
    fflush(pdf_[fld]);
  }
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
  static time_t tStart,tFinish; // Keep track of elapsed clock time

  static int first=1;
  if(first) // At beginning of run output run parameters
  {
    /* info_0.dat */
    sprintf(name_,"info%s",ext_);
    info_=fopen(name_,mode_);

    fprintf(info_,"--------------------------\n");
    fprintf(info_,"Model Specific Information\n");
    fprintf(info_,"--------------------------\n");
    modelinfo(info_);

    fprintf(info_,"\n--------------------------\n");
    fprintf(info_,"General Program Information\n");
    fprintf(info_,"-----------------------------\n");
    fprintf(info_,"Grid size=%d^%d\n",N,NDIMS);
    fprintf(info_,"Number of fields=%d\n",nflds);
    fprintf(info_,"L=%f\n",L);
    fprintf(info_,"dt=%f, dt/dx=%f\n",dt,dt/dx);
    fprintf(info_,"\n--------------------------\n");
    fprintf(info_,"Detailed Program Setup\n");
    fprintf(info_,"-----------------------------\n");
    if(hankel)
      fprintf(info_,"Hankel Fluctuation Initial Condition\n");
    else
      fprintf(info_,"EXP Fluctuation Initial Condition\n");
    fprintf(info_,"smeansvars=%d\n",smeansvars);
    fprintf(info_,"sexpansion=%d\n",sexpansion);
    fprintf(info_,"senergy=%d\n",senergy);
    fprintf(info_,"sspectra=%d\n",sspectra);
    fprintf(info_,"spdf=%d\n",spdf);
    if(spdf)
      fprintf(info_,"pdfrange=%d\n",pdfrange);

    time(&tStart);
    fprintf(info_,"\nRun began at %s",ctime(&tStart));

    /* python.dat */
    python_=fopen("python.dat","w");
    fprintf(python_,"%d %d %f %e %d %d %d %d %d %d\n",MODEL,N,L,m,smeansvars,sexpansion,senergy,sspectra,spdf,pdfrange);
    fprintf(python_,"%f %f ",initfield[0],ad);

    first=0;
  }
  else // If not at beginning record elapsed time for run
  {
    /* info_0.dat */
    time(&tFinish);
    fprintf(info_,"Run ended at %s",ctime(&tFinish));
    fprintf(info_,"\nRun from t=%f to t=%f took ",t0,t);
    readable_time((int)(tFinish-tStart),info_);
    fprintf(info_,"\n");

    /* python.dat */
    fprintf(python_,"%f %f\n",a,ad);
  }

  fflush(info_);
  fflush(python_);
  return;
}

// Calculate and save quantities (means, variances, etc.). If force>0 all infrequent calculations will be performed
void save(int force, ARG)
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


  /* Frequent calculation */
  if(smeansvars)
    meansvars(flush,flist);
  if(sexpansion)
    scale(flush);
  if(senergy)
    energy(flist);

  /* Infrequent calculation -- spectra */
  if(sspectra && flush && t>=tspectra)
  {
    spectra(flist);
    if(spdf)
      pdf(flist);
  }

  /*
  if(scheckpoint && flush && t>=tcheckpoint)
    checkpoint();
  */
}





