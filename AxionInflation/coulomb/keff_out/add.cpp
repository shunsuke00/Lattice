/*
This C++ file is calculating keff or keff/hubble_ini
in the order of keff in "initialize.cpp".
*/

#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
#include "add.hpp"

// output value is 0 -> keff, 1 -> keff/hubble=x0
#define OUTPUT 1

const double pi = (double)(2.*asin(1.));
inline double pw2(double x) {return x*x;}


int main() {
  FILE *keff_;
  char name_[100];
  char mode_[10];
  #if OUTPUT == 0
  sprintf(name_,"keff_list.dat");
  #elif OUTPUT == 1
  sprintf(name_,"x0_list.dat");
  #endif
  sprintf(mode_,"w");
  keff_=fopen(name_,mode_);

#if OUTPUT == 1
  //double phi = 5.0;
  //double phid = 0.80655;
  //double hubble = sqrt((pw2(phi)+pw2(phid))/6.);

  const double couf = 1./3.3;
  const double coef = 0.996;
  
  /*double phi = 4.6205366;
  double phid = 0.092055699;*/
  /*double phi = 4.5929167;
  double phid = 0.12242202;*/
  /*double phi = 4.5042820;
  double phid = 0.24034564;*/
  double phi = 4.4273861;
  double phid = 0.36180591;
  /*double phi = 4.3387189;
  double phid = 0.51552635;*/

  double hubble = sqrt((pw2(phid)+pow(phi,2.)+2.*coef*couf*phi*sin(phi/couf))/6.);
  printf("H0=%f\n",hubble);
#endif
  
  int px,py,pz;
  int i,j,k;
  double keff,x0;
  int N = 64;
  double L = 3.5;
  double dx = L/(double)N;
  for(i=0;i<N;i++)
    {
      px=(i<=N/2 ? i : i-N); // x-component of momentum of modes at x=i
      for(j=0;j<N;j++)
      {
        py=(j<=N/2 ? j : j-N); // y-component of momentum of modes at y=j
        // Set all modes with 0<k<N/2. The complex conjugates of these modes do not need to be set.
        for(k=1;k<N/2;k++)
        {
          pz=k; // z-component of momentum of modes at z=k
          keff=sqrt(pw2(sin(2.*pi*(double)px/(double)N))+pw2(sin(2.*pi*(double)py/(double)N))+pw2(sin(2.*pi*(double)pz/(double)N)))/dx;
          #if OUTPUT == 0
          fprintf(keff_,"%14.12f\n",keff);
          #elif OUTPUT == 1
          x0 = keff/hubble;
          fprintf(keff_,"%14.12f\n",x0);
          #endif
        }

        // Set modes with k=0 or N/2.
        if(j>N/2 || (i>N/2 && (j==0 || j==N/2))) // The complex conjugates of these modes appear explicitly on the lattice and must be set to satisfy f(-p)=f(p)*
        {
          // k=0
          keff=sqrt(pw2(sin(2.*pi*(double)px/(double)N))+pw2(sin(2.*pi*(double)py/(double)N)))/dx;
          #if OUTPUT == 0
          fprintf(keff_,"%14.12f\n",keff);
          #elif OUTPUT == 1
          fprintf(keff_,"%14.12f\n",keff/hubble);
          #endif
          // k=N/2
          keff=sqrt(pw2(sin(2.*pi*(double)px/(double)N))+pw2(sin(2.*pi*(double)py/(double)N)))/dx;
          #if OUTPUT == 0
          fprintf(keff_,"%14.12f\n",keff);
          #elif OUTPUT == 1
          fprintf(keff_,"%14.12f\n",keff/hubble);
          #endif
        }
        else if((i==0 || i==N/2) && (j==0 || j==N/2)) // The 8 "corners" of the lattice are set to real values
        {
          // k=0
          keff=0.;
          #if OUTPUT == 0
          fprintf(keff_,"%14.12f\n",keff);
          #elif OUTPUT == 1
          fprintf(keff_,"%14.12f\n",keff/hubble);
          #endif
          // k=N/2
          keff=0.;
          #if OUTPUT == 0
          fprintf(keff_,"%14.12f\n",keff);
          #elif OUTPUT == 1
          fprintf(keff_,"%14.12f\n",keff/hubble);
          #endif
        }

        fflush(keff_);
      } // End of loop over j (y-index on lattice)
    } // End of loop over i (x-index on lattice)
}
