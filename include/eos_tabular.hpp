#ifndef EOS_TABULAR_HPP
#define EOS_TABULAR_HPP

#include "eos.hpp"
#include <string>
#include <fstream>
#include <vector> 
#include <filesystem>
#include "../libInterpolate/Interpolate.hpp"
#include "../boost/function.hpp"
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>

class EOS_Tabular : public EOS {
    public:
        EOS_Tabular(std::string);
        //Copy constructor
        EOS_Tabular(const EOS_Tabular&);
        //Move constructor
        //EOS_Tabular(EOS_Tabular&&);
        //EOS_Tabular();
        EOS_Tabular(EOS_Tabular&&);
        ~EOS_Tabular();

    public:
        double p_surface;
        double e_min;
        double e_max;
        double p_max;
        double energy_from_pressure(const double&);
        double pc_from_ec(const double&);
        double dedp(const double&);
        int tab_size;

    private:
        void remove_leading_zero(std::vector<double>&);
        const std::filesystem::path eos_path = "eos_tables/";
        const std::string eos_suffix = ".dat";
        const std::string eos_prefix = "LALSimNeutronStarEOS_";
        double p, e; //tabulated pressure and energy
        const double rel_dp = 1e-5;
        // std::vector<double> log_e_tab, log_p_tab;
        // boost::math::interpolators::pchip<std::vector<double>>* e_p;
        // boost::math::interpolators::pchip<std::vector<double>>* p_e;
        //boost::function<double (double)> e_p;
        //boost::function<double (double)> p_e;     
        //std::function<double (double)> e_p;
        //std::function<double (double)> p_e;
        void initialize_splines(double*, double*);
        gsl_interp_accel* acc_e;
        gsl_interp_accel* acc_p;
        gsl_spline* e_of_p;
        gsl_spline* p_of_e;
        double e_p(double);
        double p_e(double);

        
        //_1D::CubicSplineInterpolator<double> e_p;
        //_1D::CubicSplineInterpolator<double> p_e;
};



#endif