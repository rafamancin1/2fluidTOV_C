#include "../include/eos_tabular.hpp"
#include <iostream>
#include <cmath>
#include "../include/spline.h"
#include "../boost/math/interpolators/pchip.hpp"


using boost::math::interpolators::pchip; //Piecewise Cubic Hermite interpolation

//Copy constructor 
EOS_Tabular::EOS_Tabular(const EOS_Tabular& other) : EOS_Tabular{other.eos_name} {};
// Project must always be directly under home
EOS_Tabular::EOS_Tabular(std::string eos_name) : EOS(eos_name, "Tabular") {
    std::string filename = eos_path / (eos_prefix + eos_name + eos_suffix);
    std::cout << filename << std::endl;
    std::filesystem::path home = getenv("HOME");
    std::filesystem::path project_path = home / "2fluidTOV_C";
    std::filesystem::path full_filename = project_path / filename;
    std::ifstream f_eos;
    //std::cout << "Opening file " << filename << std::endl;
    f_eos.open(full_filename);
    if (!f_eos.is_open()) {
        std::cout << "Cannot open file " << full_filename << std::endl;
        exit(0);
    };
    f_eos.seekg(0);
    std::vector<double> log_e_tab;
    std::vector<double> log_p_tab;
    std::vector<double> log_e_tab_copy;
    std::vector<double> log_p_tab_copy;
    //std::cout << "Reading file" << std::endl;
    while(true) {
        f_eos >> p;
        f_eos >> e;
        if (f_eos.eof()) { break; }
        log_e_tab.push_back(log10(e));
        log_p_tab.push_back(log10(p));
        log_e_tab_copy.push_back(log10(e));
        log_p_tab_copy.push_back(log10(p));
    }
    //remove_leading_zero(log_e_tab);
    //remove_leading_zero(log_p_tab);
    p_surface = pow(10, log_p_tab[0]);
    e_min = pow(10, log_e_tab[0]);
    e_max = pow(10, log_e_tab.back());
    p_max = pow(10, log_p_tab.back());
    f_eos.close();
    //std::vector<double> log_e_tab_copy = log_e_tab;
    //std::vector<double> log_p_tab_copy = log_p_tab;
    //e_p = pchip(std::move(log_p_tab_copy), std::move(log_e_tab_copy));
    //p_e = pchip(std::move(log_e_tab), std::move(log_p_tab));
    //tk::spline e_of_p(log_p_tab, log_e_tab);
    //tk::spline p_of_e(log_e_tab, log_p_tab);
    e_p = tk::spline(log_p_tab, log_e_tab); //e_of_p;
    p_e = tk::spline(log_e_tab, log_p_tab); //p_of_e;
    //e_p.setData(log_p_tab, log_e_tab);
    //p_e.setData(log_e_tab, log_p_tab);
};

void EOS_Tabular::remove_leading_zero(std::vector<double>& x_tab) {
    if (x_tab[0] == 0) {
        x_tab.erase(x_tab.begin());
    }
}

double EOS_Tabular::energy_from_pressure(const double& pressure) {
    // auto spline = pchip<decltype(log_p_tab)>(log_p_tab, log_e_tab);
    double e_tab = pressure < p_surface ? 0.0 : pow(10.0, e_p(log10(pressure)));
    return e_tab; 
};

double EOS_Tabular::pc_from_ec(const double& central_energy) {
    // auto spline = pchip<decltype(log_e_tab)>(log_e_tab, log_p_tab);
    const double log_energy = log10(central_energy);
    const double exponent = p_e(log_energy);
    double p_tab = pow(10.0, exponent);
    return p_tab;
};

double EOS_Tabular::dedp(const double& pressure) {
    
    if (pressure == 0) {
        return 0;
    }
    else {
        double dp = pressure * rel_dp;
        double p_upper = pressure + dp;
        double p_lower = pressure - dp;

        double eps_upper = energy_from_pressure(p_upper);
        double eps_lower = energy_from_pressure(p_lower);

        double dedp_value = (eps_upper - eps_lower) / (2 * dp);
        return dedp_value;
    }
    
   //return e_p.derivative(log10(pressure));
};