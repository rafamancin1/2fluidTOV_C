#ifndef TOV_FAMILY_HPP
#define TOV_FAMILY_HPP

#include "twofluid_TOV.hpp"
#include <boost/function.hpp>
#include "../libInterpolate/Interpolate.hpp"

template<typename T>
std::vector<double> linspace(T, T, int);
const double K_MIN = 1e6;
const double K_MAX = 1e8;
const double F_CHI_MIN = 0.01;
const double F_CHI_MAX = 0.5;
const double M_MIN = 0.2;
const double M_MAX = 2.4;
const int N_SAMPLE = 200;

class TOV_Family {
    public:
        TOV_Family(EOS_Tabular&, EOS_Poly&, const int);
        TOV_Family(EOS_Tabular&, EOS_Poly&, std::string filename);
        TOV_Family(EOS_Tabular&, EOS_Poly&, std::vector<double>&, std::vector<double>&);
        TOV_Family(EOS_Tabular&, EOS_Poly&); // For data generation
        TOV_Family();
        //TOV_Family(std::string); // For python wrapping
        double radius_from_mass(double);
        double k2_from_mass(double);
        double lambda_from_mass(double);
        double calc_lambda(const double, const double, const double);
        double calc_lambda_normalized(double, double, double);
        double calc_lambda_parallel(const double, const double, const double, const int);
        void write_to_file();
        void generate_lambda_f_points();
        std::vector<double> Rs, Ms, k2s;
        EOS_Tabular& eos1;
        EOS_Poly& eos2;   
    private:
        //boost::function<double (double)> r_m;
        //boost::function<double (double)> k2_m;
        _1D::CubicSplineInterpolator<double> r_m;
        _1D::CubicSplineInterpolator<double> k2_m;
        std::vector<double> e1s;
        std::vector<double> e2s;
        void sort_mass();
        void remove_equal_entries();
        void initialize_splines(const std::vector<double>&, const std::vector<double>&);
        double calc_e2_from_F_chi(const double&, const double&, const double&);
        double calc_F_chi(const double&, const double&, const double&);
        void reset_state();
        void add_to_vector(double, double, int, int);

        static double normalize(const double X, const double X_min, const double X_max) {return (X-X_min)/(X_max-X_min);};
        static double normalize_inverse(const double X, const double X_min, const double X_max) {return (X_max-X_min)*X + X_min;};

        int n_samples;
};

#endif