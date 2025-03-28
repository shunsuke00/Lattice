#include "parameter.hpp"
#include "physical_quantity.hpp"
#include "evolution.hpp"


int main() {
    std::ofstream _status("../status.dat", std::ios::trunc);
    if( !_status ){ std::cerr << "Failed to open 'status.txt'" << std::endl; exit(1); }
    else
    {
        _status << "# t_pr   N      a     H/m         f              df/m        ";
        _status << "energy_pr         E        epsilon       eta";
        #if INT==2
        _status << "        xi       condition1";
        #endif
        _status << std::endl;
        _status << "# -----------------------------------------------------------";
        _status << "------------------------------------------------------";
        #if INT==2
        _status << "---------------------";
        #endif
    }

    _status << std::endl;

    _status << std::noshowpos << std::fixed <<std::setprecision(3);
        _status << t << "  " << N(a);
        _status << "  " << a ;
        _status << std::noshowpos << std::fixed <<std::setprecision(5);
        _status << "  " << pow(a, s-A-1)*da;

        _status << std::showpos << std::scientific << std::setprecision(7);
        _status << "  " << f << "  " << pow(a, s-A)*df;
        _status << std::showpos << std::scientific << std::setprecision(3);
        _status << "  " << energy(f,df,a,da) << "  " << energyCons(f,df,a,da);
        _status << "  " << epsilon(f,df,a,da) << "  " << eta(f,df,a,da);
    
    #if INT==2
    _status << "  " << xi(df,a,da) << "  " << BRcondition1(df,a,da);
    #endif

    switch(METHOD)
    {
        case 1: std::cout << "Integration method is RK4" << std::endl; break;
        case 2: std::cout << "Integration method is Euler" << std::endl; break;
    }
    std::cout << "START COMPUTATION" << std::endl;

    int max_loop = total_step/output_step;

    for( int loop = 1; loop <=max_loop; ++loop )
    {
        run(output_step);
        t = loop*(dt*output_step);
        

        _status << std::endl;

        _status << std::noshowpos << std::fixed <<std::setprecision(3);
        _status << t << "  " << N(a);
        _status << "  " << a ;
        _status << std::noshowpos << std::fixed <<std::setprecision(5);
        _status << "  " << pow(a, s-A-1)*da;

        _status << std::showpos << std::scientific << std::setprecision(7);
        _status << "  " << f << "  " << pow(a, s-A)*df;
        _status << std::showpos << std::scientific << std::setprecision(3);
        _status << "  " << energy(f,df,a,da) << "  " << energyCons(f,df,a,da);
        _status << "  " << epsilon(f,df,a,da) << "  " << eta(f,df,a,da);

        #if INT==2
        _status << "  " << xi(df,a,da) << "  " << BRcondition1(df,a,da);
        #endif
    }

    std::cout << "FINISH COMPUTATION" << std::endl;
}
