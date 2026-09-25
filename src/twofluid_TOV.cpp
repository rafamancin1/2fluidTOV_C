#include "../include/twofluid_TOV.hpp"
#include "../include/conversions.hpp"
#include "../include/gsl_function_pp.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <string>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_roots.h>

#include "../include/eos_analytic.hpp"

namespace {
// Turns GSL's abort-on-error handler off while alive, so GSL failures come back
// as status codes (and become NaN results) instead of killing the process
struct GslErrorHandlerOff {
    gsl_error_handler_t* previous;
    GslErrorHandlerOff() : previous(gsl_set_error_handler_off()) {}
    ~GslErrorHandlerOff() { gsl_set_error_handler(previous); }
};

bool all_finite(const tov_state& y) {
    return std::all_of(y.begin(), y.end(), [](double v) { return std::isfinite(v); });
}
}

TwoFluid_TOV::TwoFluid_TOV(EOS_Tabular& eos1, EOS_Analytic& eos2)
    : eos1(eos1), eos2(eos2), M_B(0), R_B(0), M_D(0), R_D(0), M(0), R(0),
      y_r(0), k2(0), lambda(0), F_chi(0), dr_min(1e-2),
      ps1(0), ps2(0), active1(false), active2(false) {};

TOV_result TwoFluid_TOV::nan_result() {
    const double nan = std::numeric_limits<double>::quiet_NaN();
    return TOV_result{nan, nan, nan, nan, nan, nan, nan, nan, nan};
}

int TwoFluid_TOV::twofluid_tov_eqns(double const r, const double* const y, double* const dydr) {
    const double m1 = y[0];
    const double m2 = y[1];
    // A fluid past its surface (or an RK trial stage below it) contributes no pressure or energy
    const double p1 = (active1 && y[2] > ps1) ? y[2] : 0.0;
    const double p2 = (active2 && y[3] > ps2) ? y[3] : 0.0;
    const double y_r = y[4];

    const double e1 = eos1.energy_from_pressure(p1);
    const double e2 = eos2.energy_from_pressure(p2);
    const double dedp1 = eos1.dedp(p1);
    const double dedp2 = eos2.dedp(p2);
    const double e1_p1 = e1 + p1;
    const double e2_p2 = e2 + p2;

    const double ep1 = dedp1 * e1_p1;
    const double ep2 = dedp2 * e2_p2;
    const double p_total = p1 + p2;
    const double e_total = e1 + e2;
    const double m = m1 + m2;

    const double sdenom = r - 2 * m;
    const double denom = r * sdenom;
    const double four_pi = 4 * M_PI;
    const double r_squared = pow(r, 2);
    const double r_cubed = pow(r, 3);
    const double four_pi_r_squared = four_pi*r_squared;
    const double four_pi_r_cubed = four_pi*r_cubed;
    const double num = (four_pi_r_cubed * p_total + m);

    const double Q2nd = 4 * pow((m + four_pi_r_cubed * p_total) / denom, 2);
    const double Q1st = four_pi * r * (5 * e_total + 9 * p_total + ep1 + ep2 - 6 / (four_pi_r_squared)) / sdenom;
    const double Qr = Q1st - Q2nd;
    const double Fr = (r - four_pi_r_cubed * (e_total - p_total)) / sdenom;

    dydr[0] = four_pi_r_squared * e1;
    dydr[1] = four_pi_r_squared * e2;
    dydr[2] = -(e1_p1) * num / denom;
    dydr[3] = -(e2_p2) * num / denom;
    dydr[4] = -(pow(y_r, 2) + y_r * Fr + r_squared * Qr) / r;

    return GSL_SUCCESS;
};

// Step length h in (0, h_acc] from (r_prev, y_prev) at which the pressure of `fluid`
// (0 = baryons, 1 = DM) reaches its surface value; NaN if it cannot be located
double TwoFluid_TOV::locate_surface(int fluid, double r_prev, const tov_state& y_prev, double h_acc,
                                    gsl_odeiv2_step* stepper, const gsl_odeiv2_system* sys) {
    const int idx = 2 + fluid;
    const double p_s = fluid == 0 ? ps1 : ps2;
    tov_state y_trial(y_prev.size()), y_err(y_prev.size());
    auto excess_pressure = [&](double h) {
        y_trial = y_prev;
        if (gsl_odeiv2_step_apply(stepper, r_prev, h, &y_trial[0], &y_err[0], NULL, NULL, sys) != GSL_SUCCESS) {
            return std::numeric_limits<double>::quiet_NaN();
        }
        return y_trial[idx] - p_s;
    };
    // evolve_apply found the pressure at or below p_s at the step end. When the step barely
    // crossed, recomputing it with h_acc = r - r_prev (which carries rounding error) can land
    // just above p_s; the crossing then lies within rounding of the step end, so take it there.
    const double excess_at_end = excess_pressure(h_acc);
    if (std::isnan(excess_at_end)) { return excess_at_end; }
    if (excess_at_end >= 0.0) { return h_acc; }
    gsl_function_pp<decltype(excess_pressure)> Fp(excess_pressure);
    gsl_function* F = static_cast<gsl_function*>(&Fp);
    gsl_root_fsolver* s = gsl_root_fsolver_alloc(gsl_root_fsolver_brent);
    double h = std::numeric_limits<double>::quiet_NaN();
    if (gsl_root_fsolver_set(s, F, 0.0, h_acc) == GSL_SUCCESS) {
        int status = GSL_CONTINUE;
        for (int iter = 0; iter < 100 && status == GSL_CONTINUE; iter++) {
            if (gsl_root_fsolver_iterate(s) != GSL_SUCCESS) { break; }
            status = gsl_root_test_interval(gsl_root_fsolver_x_lower(s), gsl_root_fsolver_x_upper(s), 1e-9, 1e-12);
        }
        // The upper end has the pressure at or just below the surface value
        if (status == GSL_SUCCESS) { h = gsl_root_fsolver_x_upper(s); }
    }
    gsl_root_fsolver_free(s);
    return h;
}

TOV_result TwoFluid_TOV::integrate_two_fluid_tov(double e01, double e02) {
    GslErrorHandlerOff gsl_guard;
    reset_state();
    auto fail = [&](const std::string& why) {
        std::cerr << "[WARN] TwoFluid_TOV: " << why << " (e1=" << e01 << ", e2=" << e02
                  << " GeV/fm^3); returning NaN" << std::endl;
        reset_state();
        return nan_result();
    };
    // Central energy densities come in GeV/fm^3
    if (!std::isfinite(e01) || !std::isfinite(e02) || e01 <= 0.0 || e02 < 0.0) {
        return fail("invalid central energy densities");
    }
    const double e1c = e01 * CONVERSION::dens_GeV_fm3_to_geom;
    const double e2c = e02 * CONVERSION::dens_GeV_fm3_to_geom;

    //Initialize variables
    double dr = 0.1;
    const double p01 = eos1.pc_from_ec(e1c); // for tabular EOSs this also sets eos1.p_surface
    const double p02 = eos2.pc_from_ec(e2c);
    if (!std::isfinite(p01) || !std::isfinite(p02)) {
        return fail("central pressure is not finite (outside the EOS table?)");
    }
    ps1 = eos1.p_surface;
    ps2 = eos2.p_surface;
    const double four_pi = 4*M_PI;
    // Second order Taylor expansion about the centre: p_i(r) = p_i0 - (2pi/3)(e_i0+p_i0)(e0+3p0) r^2
    const double common_factor_p = (-four_pi/6)*(3*(p01+p02)+(e1c+e2c))*pow(dr,2);
    const double common_factor_m = (four_pi/3)*pow(dr,3);
    tov_state y = {
        e1c * common_factor_m,
        e2c * common_factor_m,
        p01 + (e1c + p01) * common_factor_p,
        p02 + (e2c + p02) * common_factor_p,
        y_0r
    };
    // A fluid whose central pressure is already below its surface value is absent (e.g. e02 = 0)
    active1 = y[2] > ps1;
    active2 = y[3] > ps2;
    if (!active1) { y[0] = 0.0; y[2] = 0.0; }
    if (!active2) { y[1] = 0.0; y[3] = 0.0; }
    double r = dr;
    // Integrator
    const int sys_dim = 5;
    gsl_odeiv2_system sys = {twofluid_tov_eqns_gsl, NULL, sys_dim, reinterpret_cast<void *>(::std::addressof(*this))};
    const gsl_odeiv2_step_type* stepper_type = gsl_odeiv2_step_rkf45;
    gsl_odeiv2_step* stepper = gsl_odeiv2_step_alloc(stepper_type, sys_dim);
    gsl_odeiv2_control* stepper_control = gsl_odeiv2_control_y_new(0.0, 1e-10); // relative error 1e-10
    gsl_odeiv2_evolve* ode_ev = gsl_odeiv2_evolve_alloc(sys_dim);
    const int max_steps = 1000000;
    const double r_max = 1e10; // m; only a safety net, physical surfaces lie far inside
    std::string failure;
    int i = 0;
    while (active1 || active2) {
        if (i++ >= max_steps) { failure = "maximum number of steps reached"; break; }
        if (r >= r_max) { failure = "surface not reached before r_max"; break; }
        const double r_prev = r;
        const tov_state y_prev = y;
        const int status = gsl_odeiv2_evolve_apply(ode_ev, stepper_control, stepper, &sys, &r, r_max, &dr, &y[0]);
        if (status != GSL_SUCCESS || !all_finite(y) || y[0] < 0 || y[1] < 0) {
            failure = "integration failed or produced unphysical values at r=" + std::to_string(r_prev) + " m";
            break;
        }
        // Did a surface lie inside the accepted step [r_prev, r]?
        const bool crossed1 = active1 && y[2] <= ps1;
        const bool crossed2 = active2 && y[3] <= ps2;
        if (!crossed1 && !crossed2) { continue; }
        const double h_acc = r - r_prev;
        const double inf = std::numeric_limits<double>::infinity();
        const double h1 = crossed1 ? locate_surface(0, r_prev, y_prev, h_acc, stepper, &sys) : inf;
        const double h2 = crossed2 ? locate_surface(1, r_prev, y_prev, h_acc, stepper, &sys) : inf;
        if (std::isnan(h1) || std::isnan(h2)) { failure = "could not locate a fluid surface"; break; }
        // Restart from the first surface; the other fluid is picked up on a later step
        const double h = std::min(h1, h2);
        tov_state y_err(sys_dim);
        y = y_prev;
        if (gsl_odeiv2_step_apply(stepper, r_prev, h, &y[0], &y_err[0], NULL, NULL, &sys) != GSL_SUCCESS || !all_finite(y)) {
            failure = "step to a fluid surface failed";
            break;
        }
        r = r_prev + h;
        if (h1 <= h2) { M_B = y[0]; R_B = r; active1 = false; y[2] = 0.0; }
        else          { M_D = y[1]; R_D = r; active2 = false; y[3] = 0.0; }
        gsl_odeiv2_evolve_reset(ode_ev);
        gsl_odeiv2_step_reset(stepper);
    };
    gsl_odeiv2_evolve_free(ode_ev);
    gsl_odeiv2_control_free(stepper_control);
    gsl_odeiv2_step_free(stepper);
    if (!failure.empty()) { return fail(failure); }

    y_r = y[4]; // at the outer surface
    M = M_B + M_D;
    R = std::max(R_D,R_B);
    F_chi = M_D / M;
    calc_k2();
    calc_lambda();
    if (!std::isfinite(F_chi) || !std::isfinite(k2) || !std::isfinite(lambda)) {
        return fail("invalid compactness, k2 or Lambda (M=" + std::to_string(M) + " m, R=" + std::to_string(R) + " m)");
    }
    TOV_result result = {.M_B=M_B, .M_D=M_D,
                         .R_B=R_B, .R_D=R_D,
                         .M=M, .R=R,
                         .F_chi=F_chi, .k2=k2,
                         .lambda=lambda};
    return result;
};

// k2 is NaN when it cannot be computed (non-positive M or R, C >= 0.5, bad denominator)
void TwoFluid_TOV::calc_k2() {
    k2 = std::numeric_limits<double>::quiet_NaN();
    if (!(M > 0 && R > 0)) { return; }
    const double C = M / R;
    if (C >= 0.5) { return; }
    const double num = ((8.0 / 5.0) * pow(1.0 - 2.0 * C, 2.0) *
                        pow(C, 5.0) * (2.0 * C * (y_r - 1.0) - y_r + 2.0));
    const double denom = (2.0 * C * (4.0 * (y_r + 1.0) * pow(C, 4.0) + (6.0 * y_r - 4.0) * pow(C, 3.0) +
             (26.0 - 22.0 * y_r) * pow(C, 2.0) + 3.0 * (5.0 * y_r - 8.0) *
             C - 3.0 * y_r + 6.0) - 3.0 * pow(1.0 - 2.0 * C, 2.0) *
             (2.0 * C * (y_r - 1.0) - y_r + 2.0) *
             log(1.0 / (1.0 - 2.0 * C)));
    if (denom > 0) { k2 = num / denom; }
}

// Lambda is NaN when k2 is invalid or the result would be negative
void TwoFluid_TOV::calc_lambda() {
    lambda = std::numeric_limits<double>::quiet_NaN();
    if (!(k2 > 0 && M > 0 && R > 0)) { return; }
    const double C = M / R;
    const double value = (2.0 / 3.0) * k2 / pow(C, 5);
    if (value >= 0) { lambda = value; }
}

void TwoFluid_TOV::reset_state() {
    M_B = 0;
    M_D = 0;
    R_B = 0;
    R_D = 0;
    M = 0;
    R = 0;
    y_r = 0;
    k2 = 0;
    lambda = 0;
    F_chi = 0;
    dr_min = 1e-2;
    active1 = false;
    active2 = false;
};

void TwoFluid_TOV::print_result(TOV_result& result) {
    const double M_B_Msun = result.M_B * CONVERSION::mass_geom_to_Msun;
    const double M_D_Msun = result.M_D * CONVERSION::mass_geom_to_Msun;
    const double R_B_km = result.R_B * CONVERSION::m_to_km;
    const double R_D_km = result.R_D * CONVERSION::m_to_km;
    const double M_Msun = result.M * CONVERSION::mass_geom_to_Msun;
    const double R = result.R * CONVERSION::m_to_km;
    const double lambda = result.lambda;

    std::cout << "Baryon mass(M_sun) = " << M_B_Msun << std::endl;
    std::cout << "DM mass(M_sun) = " << M_D_Msun << std::endl;
    std::cout << "Baryon radius(km) = " << R_B_km << std::endl;
    std::cout << "DM radius(km) = " << R_D_km << std::endl;
    std::cout << "Total mass(M_sun) = " << M_Msun << std::endl;
    std::cout << "Radius(km) = " << R << std::endl;
    std::cout << "Lambda = " << lambda << std::endl;
    std::cout << "F_chi = " << result.F_chi << std::endl;
};
