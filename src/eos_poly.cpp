#include "../include/eos_poly.hpp"
#include <boost/bind.hpp>
#include <functional>
#include <boost/math/tools/minima.hpp>

using boost::math::tools::brent_find_minima; //Brent's method to locate minima of function 

EOS_Poly::EOS_Poly(std::string eos_name, const double K, const double Gamma) : EOS(eos_name, "Polytropic"), K(K), Gamma(Gamma), p_surface(1e-35) {};
EOS_Poly::EOS_Poly(const EOS_Poly& other) : EOS_Poly{other.eos_name, other.K, other.Gamma} {};

double EOS_Poly::pc_from_ec(const double& central_energy) {
    /*this function is called only to find the central pressure given the value of central energy density  
    thus we must use root-finding to find the actual value since central energy density epsilon is different from central density rho*/
    double p_c_approx = K*pow(central_energy,Gamma); //set approximate value for central pressure
    const int double_bits = std::numeric_limits<double>::digits; // set numerical accuracy of double type for Brent's method
    /* using Brent's method alone is not enough to locate the desired minima, since the function pressure_minima has in fact two roots.
    For that reason we must find the first root of brent's method and to do that we check whether the function at the root obtained 
    is decreasing, in which case we lower the upper bound (upperCoef) of the method's bracketing to gradually reach the first root*/
    double upperCoef = 15.0; 
    double lowerCoef = 0.1;
    p_c = brent_find_minima(boost::bind(&EOS_Poly::pressure_minima,this,_1, central_energy), lowerCoef*p_c_approx,upperCoef*p_c_approx,double_bits);
    double pressureCheck = p_c.first+0.005;
    bool BrentWentBad = (pressureCheck-K*pow(central_energy-pressureCheck/(Gamma-1.0),Gamma))<0 || fabs(p_c.second) || std::isnan(p_c.second);
    while(BrentWentBad && upperCoef>=0.1){
      upperCoef-=0.01;
      p_c = brent_find_minima(boost::bind(&EOS_Poly::pressure_minima,this,_1, central_energy), lowerCoef*p_c_approx,upperCoef*p_c_approx,double_bits);
      pressureCheck = p_c.first+0.005;
      BrentWentBad = (pressureCheck-K*pow(central_energy-pressureCheck/(Gamma-1.0),Gamma))<0 || fabs(p_c.second)>1e-5 || std::isnan(p_c.second);
    }
    return p_c.first;
}