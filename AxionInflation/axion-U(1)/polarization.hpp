#ifndef _POLARIZATION_
#define _POLARIZATION_

#include "parameters.hpp"

const double dx=L/(double)N;  //NON-ADJUSTABLE VARIABLES

/* Polarization vector (positive and K argument) */
/*--------------------------------------------------*/
#if POLARIZATION==0
inline double pola_x_re(double px, double py, double pz)
{
  return ( -px*pz/sqrt(2.*(pw2(px)+pw2(py))*(pw2(px)+pw2(py)+pw2(pz))) );
}

inline double pola_x_im(double px, double py, double pz)
{
  return ( -py/sqrt(2.*(pw2(px)+pw2(py))) );
}

inline double pola_y_re(double px, double py, double pz)
{
  return ( -py*pz/sqrt(2.*(pw2(px)+pw2(py))*(pw2(px)+pw2(py)+pw2(pz))) );
}

inline double pola_y_im(double px, double py, double pz)
{
  return ( px/sqrt(2.*(pw2(px)+pw2(py))) );
}

inline double pola_z_re(double px, double py, double pz)
{
  return ( sqrt( (pw2(px)+pw2(py))/(2.*(pw2(px)+pw2(py)+pw2(pz))) ) );
}

inline double pola_z_im(double px, double py, double pz)
{
  return 0.;
}

/* Polarization vector (positive and keff argument) */
/*--------------------------------------------------*/
#elif POLARIZATION==1
inline double pola_x_re(double px, double py, double pz)
{
  double px_eff,py_eff,pz_eff;
  px_eff = sin(2.*pi*(double)px/(double)N)/dx;
  py_eff = sin(2.*pi*(double)py/(double)N)/dx;
  pz_eff = sin(2.*pi*(double)pz/(double)N)/dx;
  return ( -px_eff*pz_eff/sqrt(2.*(pw2(px_eff)+pw2(py_eff))*(pw2(px_eff)+pw2(py_eff)+pw2(pz_eff))) );
}

inline double pola_x_im(double px, double py, double pz)
{
  double px_eff,py_eff,pz_eff;
  px_eff = sin(2.*pi*(double)px/(double)N)/dx;
  py_eff = sin(2.*pi*(double)py/(double)N)/dx;
  pz_eff = sin(2.*pi*(double)pz/(double)N)/dx;
  return ( py_eff/sqrt(2.*(pw2(px_eff)+pw2(py_eff))) );
}

inline double pola_y_re(double px, double py, double pz)
{
  double px_eff,py_eff,pz_eff;
  px_eff = sin(2.*pi*(double)px/(double)N)/dx;
  py_eff = sin(2.*pi*(double)py/(double)N)/dx;
  pz_eff = sin(2.*pi*(double)pz/(double)N)/dx;
  return ( -py_eff*pz_eff/sqrt(2.*(pw2(px_eff)+pw2(py_eff))*(pw2(px_eff)+pw2(py_eff)+pw2(pz_eff))) );
}

inline double pola_y_im(double px, double py, double pz)
{
  double px_eff,py_eff,pz_eff;
  px_eff = sin(2.*pi*(double)px/(double)N)/dx;
  py_eff = sin(2.*pi*(double)py/(double)N)/dx;
  pz_eff = sin(2.*pi*(double)pz/(double)N)/dx;
  return ( -px_eff/sqrt(2.*(pw2(px_eff)+pw2(py_eff))) );
}

inline double pola_z_re(double px, double py, double pz)
{
  double px_eff,py_eff,pz_eff;
  px_eff = sin(2.*pi*(double)px/(double)N)/dx;
  py_eff = sin(2.*pi*(double)py/(double)N)/dx;
  pz_eff = sin(2.*pi*(double)pz/(double)N)/dx;
  return ( sqrt( (pw2(px_eff)+pw2(py_eff))/(2.*(pw2(px_eff)+pw2(py_eff)+pw2(pz_eff))) ) );
}

inline double pola_z_im(double px, double py, double pz)
{
  return 0.;
}

#endif


#endif