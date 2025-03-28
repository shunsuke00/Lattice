#ifndef _INTERACTION_
#define _INTERACTION_

#include "model.hpp"

/* coupling */
double g = 42.;
//double g = 25.;
//double g = 3.3*6.7;

/* effective coupling */
double xi(double df, double a, double da)
{
    return (std::abs(0.5*g*df*a/da));
}

/* source term for df e.o.m */
double source(double df, double a, double da)
{
    return ((-2.4)*pow(10,-4.)*g*pow(m,2.)*pow(a,2.*s+2.*A-4.)*pow(da/a,4.)*std::exp(2*PI*xi(df,a,da))/pow(xi(df,a,da),4.));
}

/* gauge field background energy */
double gauge_energy(double df, double a, double da)
{
    return (1.4*pow(10,-4.)*pow(m,2.)*pow(a,2.*s-4.*A)*pow(da/a,4.)*std::exp(2*PI*xi(df,a,da))/pow(xi(df,a,da),3.));
}

/* negligible backreaction condition (stringer one) */
double BRcondition1(double df, double a, double da)
{
    return (m*pow(a,s-3.)*pow(da,2.)*std::exp(PI*std::abs(xi(df,a,da)))/(26.*PI*std::abs(df)*pow(std::abs(xi(df,a,da)),3./2.)));
}


#endif