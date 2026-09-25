#ifndef TOV_FAMILY_HPP
#define TOV_FAMILY_HPP

#include <functional>
#include "eos_analytic.hpp"
#include "eos_tabular.hpp"
#include "eos_poly.hpp"
#include "eos_SIDM.hpp"
#include "twofluid_TOV.hpp"
#include "../libInterpolate/Interpolate.hpp"

template<typename T>
std::vector<double> linspace(T, T, int);
const double E1_MIN = 0.15; // lower bound of the baryonic central energy density grid (GeV/fm^3)
const double E1_MAX = 2.4;  // upper bound of the baryonic central energy density grid (GeV/fm^3)
const int N_SAMPLE = 400;
const double E2_MAX_DEFAULT = 2.0; // upper bound on the DM central energy density searched for a target F_chi (GeV/fm^3)

class TOV_Family {
    public:
        TOV_Family(EOS_Tabular&, EOS_Analytic&, std::vector<double>&, std::vector<double>&);
        TOV_Family(EOS_Tabular&, EOS_Analytic&, double F_chi, double e2_max = E2_MAX_DEFAULT);
        TOV_Family(EOS_Tabular&, EOS_Analytic&);
        TOV_Family(EOS_Tabular&);
        double radius_from_mass(double);
        double k2_from_mass(double);
        double lambda_from_mass(double);
        TOV_result calc_lambda_and_mass_directly(double, double);
        TOV_result calc_lambda_and_mass_directly_v2(double, double);
        double calc_lambda(const double, const double);
        double calc_e2_from_F_chi_v2(const double&, const double&);
        void set_analytic_eos(EOS_Analytic& eos2);
        std::vector<double> Rs, Ms, k2s, lambdas, RBs, RDs;
        EOS_Tabular& eos1;
        std::reference_wrapper<EOS_Analytic> eos2; // rebindable by set_analytic_eos
        std::vector<double> e1s;
        std::vector<double> e2s;
        double e2_max = E2_MAX_DEFAULT; // targets needing a larger e2 are unreachable (NaN)
    private:
        // Splines over the stable branch (masses up to the maximum mass) only
        _1D::CubicSplineInterpolator<double> r_m;
        _1D::CubicSplineInterpolator<double> k2_m;
        double M_spline_min = 0.0;
        double M_spline_max = 0.0;
        bool splines_ready = false;

        void generate_e2s_from_F_chi(double);
        void initialize_splines(const std::vector<double>&, const std::vector<double>&);
        double calc_F_chi(const double&, const double&);
        void reset_state();
        void check_splines() const;
};

#endif
