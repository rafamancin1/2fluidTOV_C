#ifndef EOS_POLY_HPP
#define EOS_POLY_HPP

#include <string>
#include <cmath>
#include "eos.hpp"


class EOS_Poly : public EOS {
    public:
        EOS_Poly(std::string, const double, const double);
        // Copy constructor
        EOS_Poly(const EOS_Poly&);
        EOS_Poly();

    public:
        double K;
        double Gamma;
        double p_surface;

    public:
        double pc_from_ec(const double&); // central pressure from central density
        double energy_from_pressure(const double& pressure) {
            if (pressure == 0) {
                return 0;
            }
            else {
                double first_term = std::pow(pressure/K, 1.0/Gamma);
                double second_term = pressure / (Gamma - 1);
                return first_term + second_term;
            }
        };
        double dedp(const double& pressure) {
            if (pressure == 0) {
                return 0;
            }
            else {
                double first_term = (1/(K*Gamma))*std::pow(pressure/K, 1.0/Gamma - 1);
                double second_term = 1.0/(Gamma - 1);
                return first_term + second_term;
            }
        };

    private:
        double pressure_minima(const double& pressure, const double& central_energy) {
            return std::fabs(pressure - K*std::pow(central_energy - pressure/(Gamma - 1.0), Gamma));
        }; // minima function used to evaluate central pressure via root-finding

    private:
        std::pair<double, double> p_c; // returned central pressure

};

#endif