#include "../include/eos_analytic.hpp"
#include "../include/eos_poly.hpp"
#include <boost/bind.hpp>
#include <functional>
#include <boost/math/tools/minima.hpp>

using boost::math::tools::brent_find_minima;

EOS_Poly::EOS_Poly(std::string eos_name, const double K, const double Gamma)
    : EOS_Analytic("Polytropic"), K(K), Gamma(Gamma), p_surface(1e-20) {}

EOS_Poly::EOS_Poly(const EOS_Poly& other)
    : EOS_Poly{other.eos_name, other.K, other.Gamma} {}

EOS_Poly::EOS_Poly()
    : EOS_Poly("null", 0.0, 0.0) {}

double EOS_Poly::pc_from_ec(const double& central_energy) const {
    double p_c_approx = K*pow(central_energy,Gamma);
    const int double_bits = std::numeric_limits<double>::digits;
    double upperCoef = 15.0;
    double lowerCoef = 0.1;
    auto self = const_cast<EOS_Poly*>(this); // p_c is mutable for caching
    self->p_c = brent_find_minima(boost::bind(&EOS_Poly::pressure_minima, self, _1, central_energy),
                                  lowerCoef*p_c_approx, upperCoef*p_c_approx, double_bits);
    double pressureCheck = self->p_c.first+0.005;
    bool BrentWentBad = (pressureCheck-K*pow(central_energy-pressureCheck/(Gamma-1.0),Gamma))<0
                        || fabs(self->p_c.second)>1e-5
                        || std::isnan(self->p_c.second);
    while(BrentWentBad && upperCoef>=0.1){
        upperCoef-=0.01;
        self->p_c = brent_find_minima(boost::bind(&EOS_Poly::pressure_minima, self, _1, central_energy),
                                      lowerCoef*p_c_approx, upperCoef*p_c_approx, double_bits);
        pressureCheck = self->p_c.first+0.005;
        BrentWentBad = (pressureCheck-K*pow(central_energy-pressureCheck/(Gamma-1.0),Gamma))<0
                       || fabs(self->p_c.second)>1e-5
                       || std::isnan(self->p_c.second);
    }
    self->p_surface = 1e-10*self->p_c.first;
    return self->p_c.first;
}