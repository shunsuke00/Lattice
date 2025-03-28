#ifndef _PHYSICAL_QUANTITY_
#define _PHYSICAL_QUANTITY_

#include "parameter.hpp"
#include "model.hpp"

#if INT==2
#include "interaction.hpp"
#endif

//physical quantity
//-----------------------------------------------------------------------
//energy density
double energy(double f, double df, double a, double da)
{
    #if INT==1
    return pow(df, 2)/(2*pow(a,2*A)) + potential(f, a);
    #elif INT==2
    if(t<0.001)
        return pow(df, 2)/(2*pow(a,2*A)) + potential(f, a);
    else
        return pow(df, 2)/(2*pow(a,2*A)) + potential(f, a) + gauge_energy(df,a,da);
    #endif
}

//pressure
double pressure(double f, double df, double a, double da)
{
    #if INT==1
    return pow(df, 2)/(2*pow(a,2*A)) - potential(f, a);
    #elif INT==2
    if(t<0.001)
        return pow(df, 2)/(2*pow(a,2*A)) - potential(f, a);
    else
        return pow(df, 2)/(2*pow(a,2*A)) - potential(f, a) + gauge_energy(df,a,da)/3.;
    #endif
}

//the quantity (must be 1) to confirm energy conservation
double energyCons(double f, double df, double a, double da)
{
    return 3 * pow(da/a, 2) / (energy(f,df,a,da) * pow(a,2*A));
}

//slow-roll parameter "epsilon"
double epsilon(double f, double df, double a, double da)
{
    return 3.*(1.+pressure(f,df,a,da)/energy(f,df,a,da))/2.;
}

double eta(double f, double df, double a, double da)
{
    //return (-3.-pow(a,1.+2.*A)*potentialDeriv(f,a)/(da*df)); // slow-roll approximation
    return (-3.*pow(a,1.-2.*A)*df*(potential(f,a)*(3.*df*da/a+pow(a,2.*A)*potentialDeriv(f,a))-0.5*pow(df,2.)*potentialDeriv(f,a))/(da*epsilon(f,df,a,da)*pow(energy(f,df,a,da),2.))); // exact fomula(single inflation)
}

//e-folding from starting program (a=1)
double N(double a)
{
    return std::log(a);
}


//time derivative of scale factor
//----------------------------------------------------------------------
double da = sqrt(energy(f,df,a,da)/3.);


#endif