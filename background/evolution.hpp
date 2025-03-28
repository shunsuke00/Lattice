#ifndef _EVOLUTION_RK4_
#define _EVOLUTION_RK4_

#include "physical_quantity.hpp"

//evolution equation
//---------------------------------------------------------------------------
//field
double equationF(double df)
{
    return df*dt;
}

//time derivative of field
double equationDF(double f, double df, double a, double da)
{
    #if INT==1
    return (-(s-A+3)*da*df/a - pow(a,2*A)*potentialDeriv(f, a)) * dt;
    #elif INT==2
    if(t<0.001)
        return (-(s-A+3)*da*df/a - pow(a,2*A)*potentialDeriv(f, a)) * dt;
    else
        return (-(s-A+3)*da*df/a - pow(a,2*A)*potentialDeriv(f, a) + source(df,a,da)) * dt;
    #endif
}

//scale factor
double equationA(double da)
{
    return da * dt;
}

//time derivative of scale factor
double equationDA(double f, double df, double a, double da)
{
    return (-s*pow(da, 2)/a - pow(a,2*A+1)*(pow(-1,A)* energy(f,df,a,da)+3*pressure(f,df,a,da))/6) * dt;
}


//Integration
//--------------------------------------------------------------------
#if METHOD == 1
    //4次のルンゲクッタ法
    void run(int output_step)
    {
        for(int n = 0; n < output_step; ++n)
        {
            double k1_f = equationF(df);
            double k1_df = equationDF(f,df,a,da);
            double k1_a = equationA(da);
            double k1_da = equationDA(f,df,a,da);

            double k2_f = equationF(df+k1_df/2);
            double k2_df = equationDF(f+k1_f/2,df+k1_df/2,a+k1_a/2,da+k1_da/2);
            double k2_a = equationA(da+k1_da/2);
            double k2_da = equationDA(f+k1_f/2,df+k1_df/2,a+k1_a/2,da+k1_da/2);

            double k3_f = equationF(df+k2_df/2);
            double k3_df = equationDF(f+k2_f/2,df+k2_df/2,a+k2_a/2,da+k2_da/2);
            double k3_a = equationA(da+k2_da/2);
            double k3_da = equationDA(f+k2_f/2,df+k2_df/2,a+k2_a/2,da+k2_da/2);

            double k4_f = equationF(df+k3_df);
            double k4_df = equationDF(f+k3_f,df+k3_df,a+k3_a,da+k3_da);
            double k4_a = equationA(da+k3_da);
            double k4_da = equationDA(f+k3_f,df+k3_df,a+k3_a,da+k3_da);

            f  += k1_f/6 + k2_f/3 + k3_f/3 + k4_f/6;
            df += k1_df/6 + k2_df/3 + k3_df/3 + k4_df/6;
            a  += k1_a/6 + k2_a/3 + k3_a/3 + k4_a/6;
            da += k1_da/6 + k2_da/3 + k3_da/3 + k4_da/6;
        }
    }
#elif METHOD == 2
    //オイラー法
    void run(int output_step)
    {
        for(int n = 0; n < output_step; ++n)
        {
            double k_f  = equationF(df);
            double k_df = equationDF(f, df, a, da);
            double k_a  = equationA(da);
            double k_da = equationDA(f, df, a, da);

            f  += k_f;
            df += k_df;
            a  += k_a;
            da += k_da;
        }
    }
#endif


#endif