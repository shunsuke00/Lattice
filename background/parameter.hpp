#ifndef ADD_HPP
#define ADD_HPP

#include <fstream>
#include <iostream>
#include <iomanip>
#include <stdio.h>
#include <numeric>
#include <cmath>

/* circle ratio */
#define PI 3.14159265359


//integration method
//"1"(rungekutta),"2"(Euler)
#define METHOD 1

//interacting field
//"1"(single field inf),"2"(axion-U(1))
#define INT 1

/* time coordinate: "0"(physical),"1"(conformal) */
double A = 1.;

/* rescaling parameter */
double s = 1.;

/* field value */
// double f  = -14.5;
// double df = 0.8152128;
// double df = -0.8233865;
//double f  = -5.5;
//double df = 0.0;
double f  = -4.8;
double df = 0.;
//double f = -4.4273861;
//double df = 0.36180591;

/* scale factor */
double a = 1.;
// double da;

/* program time */
double t = 0;

/* program time_step */
double dt = 1.0e-4;
int output_step = 2.e+1;
//int total_step = 1.2e+4;
//int total_step = 6.e+4;
int total_step = 5.0e+5;

/*
    this time_step "dt" follow "Simulating the inflationary Universe"
    this "total_step*dt" indicate  a_end = 1000 

    in this program, time coordinate is "physical time"
*/


#endif
