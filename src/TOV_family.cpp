#include "../include/TOV_family.hpp"
#include "../include/conversions.hpp"
#include "../include/gsl_function_pp.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_min.h>
#include <gsl/gsl_roots.h>

EOS_Poly null_eos("Polytropic", 0.0, 0.0);

template<typename T>
std::vector<double> linspace(T start_in, T end_in, int num_in)
{

  std::vector<double> linspaced;

  double start = static_cast<double>(start_in);
  double end = static_cast<double>(end_in);
  double num = static_cast<double>(num_in);

  if (num == 0) { return linspaced; }
  if (num == 1)
    {
      linspaced.push_back(start);
      return linspaced;
    }

  double delta = (end - start) / (num - 1);

  for(int i=0; i < num-1; ++i)
    {
      linspaced.push_back(start + delta * i);
    }
  linspaced.push_back(end); // I want to ensure that start and end
                            // are exactly the same as the input
  return linspaced;
}

TOV_Family::TOV_Family(EOS_Tabular& eos1, EOS_Analytic& eos2, std::vector<double>& e1s, std::vector<double>& e2s)
    : eos1(eos1), eos2(eos2), e1s(e1s), e2s(e2s) {
    initialize_splines(e1s, e2s);
}

// Only sets up the e1 grid; used for direct single-star solves and calc_lambda
TOV_Family::TOV_Family(EOS_Tabular& eos1, EOS_Analytic& eos2)
    : eos1(eos1), eos2(eos2) {
    e1s = linspace(E1_MIN, E1_MAX, N_SAMPLE);
}

TOV_Family::TOV_Family(EOS_Tabular& eos1, EOS_Analytic& eos2, double F_chi, double e2_max)
    : eos1(eos1), eos2(eos2), e2_max(e2_max) {
    e1s = linspace(E1_MIN, E1_MAX, N_SAMPLE);
    generate_e2s_from_F_chi(F_chi);
    initialize_splines(e1s, e2s);
}

TOV_Family::TOV_Family(EOS_Tabular& eos1) : eos1(eos1), eos2(null_eos) {
    const double e1_min = 0.3;
    const double e1_max = 1.44;
    e1s = linspace(e1_min, e1_max, N_SAMPLE);
    e2s = std::vector<double>(N_SAMPLE, 0.0);
    initialize_splines(e1s, e2s);
}

double TOV_Family::calc_F_chi(const double& e1, const double& e2) {
    TwoFluid_TOV model(eos1, eos2);
    TOV_result res = model.integrate_two_fluid_tov(e1, e2);
    return res.F_chi;
}

// DM central energy density (GeV/fm^3) in [0, e2_max] giving the target F_chi at fixed e1;
// NaN when the target is invalid, unreachable below e2_max, or the search fails
double TOV_Family::calc_e2_from_F_chi_v2(const double& F_chi, const double& e1) {
    const double nan = std::numeric_limits<double>::quiet_NaN();
    if (F_chi == 0.0) {
        return 0.0;
    }
    if (!std::isfinite(F_chi) || F_chi < 0.0 || F_chi >= 1.0) {
        return nan;
    }
    const double lower_bound = 0.0;
    const double upper_bound = e2_max;
    int status;
    double root = nan;
    auto F_chi_minima = [&F_chi, &e1, this] (double e2) {return 100*(calc_F_chi(e1, e2)-F_chi);};
    gsl_function_pp<decltype(F_chi_minima)> Fp(F_chi_minima);
    gsl_function *F = static_cast<gsl_function*>(&Fp);
    gsl_error_handler_t *default_handler = gsl_set_error_handler_off();
    // F_chi(e2) rises from 0 at e2 = 0, so the target is reachable only if F_chi(e2_max) >= target
    const double f_upper = F_chi_minima(upper_bound);
    if (std::isfinite(f_upper) && f_upper >= 0.0) {
        gsl_root_fsolver *s = gsl_root_fsolver_alloc(gsl_root_fsolver_brent);
        status = gsl_root_fsolver_set(s, F, lower_bound, upper_bound);
        // F_chi(e2) is very steep near the core-to-halo transition (small m_chi), so the
        // tolerance on e2 must be relative: an absolute one gives F_chi errors of O(0.1)
        for (int iter = 0; iter < 200 && status == GSL_SUCCESS; iter++) {
            status = gsl_root_fsolver_iterate(s);
            if (status != GSL_SUCCESS) { break; } // e.g. a failed (NaN) trial solve
            status = gsl_root_test_interval(gsl_root_fsolver_x_lower(s), gsl_root_fsolver_x_upper(s), 1e-14, 1e-10);
            if (status == GSL_SUCCESS) {
                root = gsl_root_fsolver_root(s);
                break;
            }
            status = GSL_SUCCESS; // GSL_CONTINUE: keep iterating
        }
        gsl_root_fsolver_free(s);
    }
    gsl_set_error_handler(default_handler); // to reset the error handler
    return root;
}

void TOV_Family::generate_e2s_from_F_chi(double F_chi) {
    if (F_chi == 0.0) {
        e2s = std::vector<double>(e1s.size(), 0.0);
    }
    else {
        e2s = std::vector<double>(e1s.size());
        for (size_t i = 0; i < e1s.size(); i++) {
            e2s[i] = calc_e2_from_F_chi_v2(F_chi, e1s[i]);
        }
    }
}

TOV_result TOV_Family::calc_lambda_and_mass_directly(double e1, double F_chi) {
    TwoFluid_TOV model(eos1, eos2);
    double e2 = calc_e2_from_F_chi_v2(F_chi, e1);
    if (std::isnan(e2)) {
        return TwoFluid_TOV::nan_result();
    }
    TOV_result result = model.integrate_two_fluid_tov(e1, e2);
    return result;
}

TOV_result TOV_Family::calc_lambda_and_mass_directly_v2(double e1, double F_chi) {
    TwoFluid_TOV model(eos1, eos2);
    double e2 = F_chi * e1;
    TOV_result result = model.integrate_two_fluid_tov(e1, e2);
    return result;
}

double TOV_Family::calc_lambda(const double F_chi, const double mass) {
    generate_e2s_from_F_chi(F_chi);
    initialize_splines(e1s, e2s);
    return lambda_from_mass(mass);
}

void TOV_Family::reset_state() {
    Rs.clear();
    Ms.clear();
    k2s.clear();
    lambdas.clear();
    RBs.clear();
    RDs.clear();
    splines_ready = false;
}

void TOV_Family::initialize_splines(const std::vector<double>& e1s, const std::vector<double>& e2s) {
    reset_state();
    TwoFluid_TOV model(eos1, eos2);
    for (size_t i = 0; i < e1s.size(); i++) {
        // e2 is NaN where the target F_chi was unreachable; that point is already a failure
        TOV_result res = std::isnan(e2s[i]) ? TwoFluid_TOV::nan_result() : model.integrate_two_fluid_tov(e1s[i], e2s[i]);
        Rs.push_back(res.R);
        Ms.push_back(res.M);
        k2s.push_back(res.k2);
        RBs.push_back(res.R_B);
        RDs.push_back(res.R_D);
        lambdas.push_back(res.lambda);
        model.reset_state();
    }
    // The splines need M strictly increasing, so keep only the stable branch:
    // solutions up to the maximum mass, skipping failed (NaN) or non-increasing points.
    size_t i_max = 0;
    bool any_finite = false;
    for (size_t i = 0; i < Ms.size(); i++) {
        if (std::isfinite(Ms[i]) && (!any_finite || Ms[i] > Ms[i_max])) { i_max = i; any_finite = true; }
    }
    if (!any_finite) { return; }
    std::vector<double> Ms_stable, Rs_stable, k2s_stable;
    for (size_t i = 0; i <= i_max; i++) {
        const double M_prev = Ms_stable.empty() ? 0.0 : Ms_stable.back();
        if (Ms[i] > M_prev) {
            Ms_stable.push_back(Ms[i]);
            Rs_stable.push_back(Rs[i]);
            k2s_stable.push_back(k2s[i]);
        }
    }
    if (Ms_stable.size() < 3) { return; }
    r_m.setData(Ms_stable, Rs_stable);
    k2_m.setData(Ms_stable, k2s_stable);
    M_spline_min = Ms_stable.front();
    M_spline_max = Ms_stable.back();
    splines_ready = true;
}

void TOV_Family::check_splines() const {
    if (!splines_ready) {
        throw std::runtime_error("TOV_Family: no M(R), k2(M) splines available "
                                 "(family not integrated, or fewer than 3 stable-branch solutions)");
    }
}

// Mass in M_sun; NaN outside the stable-branch mass range of the family
double TOV_Family::radius_from_mass(double mass) {
    check_splines();
    mass *= CONVERSION::Msun_to_mass_geom;
    if (mass < M_spline_min || mass > M_spline_max) { return std::numeric_limits<double>::quiet_NaN(); }
    return r_m(mass);
}

double TOV_Family::k2_from_mass(double mass) {
    check_splines();
    mass *= CONVERSION::Msun_to_mass_geom;
    if (mass < M_spline_min || mass > M_spline_max) { return std::numeric_limits<double>::quiet_NaN(); }
    return k2_m(mass);
}

double TOV_Family::lambda_from_mass(double mass) {
    const double radius = radius_from_mass(mass);
    const double k2 = k2_from_mass(mass);
    mass *= CONVERSION::Msun_to_mass_geom;
    const double C = mass / radius;
    const double lambda = (2.0 / 3.0) * k2 / pow(C, 5);
    return lambda;
}

// Rebinds the DM EOS; results computed with the previous EOS are discarded
void TOV_Family::set_analytic_eos(EOS_Analytic& eos2) {
    this->eos2 = eos2;
    reset_state();
}
