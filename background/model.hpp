#ifndef _MODEL_
#define _MODEL_

#include "parameter.hpp"

//inflaton mass
//double m = 5.1*pow(10.,-6.);
double m = 1.399*pow(10.,-7);

//potential parameters
double couf = 1./3.3;
//double coef = 0.996;
double coef = 0.9958;


//inflaton potential
double potential(double f, double a)
{
    // return pow(f, 2)/(2*pow(a, 2*s)); //chaotic inflation
    // return pow(f, 2)*(1 + 0.01*std::tanh((f-14.35)/0.005))/(2*pow(a, 2*s));   //step potential
    return (0.5*pow(f,2.)+coef*couf*f*std::sin(f/couf))/pow(a,2.*s);  //bunpy axion inflation
}

//derivative of inflaton potential
double potentialDeriv(double f, double a)
{
    // return f/pow(a, 2*s);   //chaotic inflation
    // return f*(1 + 0.01*std::tanh((f-14.35)/0.005))/pow(a, 2*s) + pow(f, 2)*0.01/(2*0.005*pow(a, 2*s)*pow(std::cosh((f-14.35)/0.005), 2));   //step potential
    return (f+coef*(couf*std::sin(f/couf)+f*std::cos(f/couf)))/pow(a, 2*s);  //bunpy axion inflation
}


#endif