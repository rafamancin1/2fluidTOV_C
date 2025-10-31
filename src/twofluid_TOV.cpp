#include "../include/twofluid_TOV.hpp"
//#include "../boost/numeric/odeint.hpp"
#include "../include/conversions.hpp"
#include <boost/bind.hpp>
#include <cmath>
#include <gsl/gsl_odeiv2.h>
#include <numeric>

#include "../include/eos_analytic.hpp"




// TwoFluid_TOV::TwoFluid_TOV(EOS_Tabular* eos1, EOS_Tabular* eos2) : eos1(eos1), eos2(eos2){};
TwoFluid_TOV::TwoFluid_TOV(EOS_Tabular& eos1, EOS_Analytic& eos2)
    : eos1(eos1), eos2(eos2), M_B(0), M_D(0), R_B(0), R_D(0), M(0), R(0), dr_min(1e-2) {};
// TwoFluid_TOV::TwoFluid_TOV(EOS_Poly* eos1, EOS_Poly* eos2) : eos1(eos1), eos2(eos2) {};

int TwoFluid_TOV::twofluid_tov_eqns(double const r, const double* const y, double* const dydr) {
    const double m1 = y[0];
    const double m2 = y[1];
    const double p1 = y[2] > eos1.p_surface ? y[2] : 0.0;
    const double p2 = y[3] > eos2.p_surface ? y[3] : 0.0;
    const double y_r = y[4];

    if (p1 == 0.0 && M_B == 0) { M_B = m1; R_B = r; }
    if (p2 == 0.0 && M_D == 0) { M_D = m2; R_D = r; }

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

/*
void TwoFluid_TOV::tov_step(controlled_runge_kutta<runge_kutta_dopri5<tov_state>>& stepper, tov_state& y, double& r, double& dr) {
    tov_state y_old = y;
    double r_old = r;
    //double dr_old = dr;
    controlled_step_result res = stepper.try_step(boost::bind(&TwoFluid_TOV::twofluid_tov_eqns, this, _1, _2, _3), y, r, dr);
    while (y[2] < 0 || y[3] < 0 || res == fail) {
        //std::cout << "Pressure went negative. Decreasing step-size" << std::endl;
        y = y_old;
        r = r_old;
        if (dr - 10 > 0) {
            dr -= 10;
        }
        else {
            dr *= 1e-1;
        }
        res = stepper.try_step(boost::bind(&TwoFluid_TOV::twofluid_tov_eqns, this, _1, _2, _3), y, r, dr);
        if (y[2] <= eos1.p_surface && y[2] != 0) {y[2] = 0;};
        if (y[3] <= eos2.p_surface && y[3] != 0) {y[3] = 0;};
        // just in case this is a possibility, which I don't know  
        // if (res == success) { res = fail;} 
    };
    
}
*/

TOV_result TwoFluid_TOV::integrate_two_fluid_tov(double e01, double e02) {
    // Some conversions first...
    // e0 in GeV/fm3
    e01 *= CONVERSION::dens_GeV_fm3_to_geom;
    e02 *= CONVERSION::dens_GeV_fm3_to_geom;


    //Initialize variables
    double dr = 0.1;
    const double p01 = eos1.pc_from_ec(e01);
    const double p02 = eos2.pc_from_ec(e02);
    const double four_pi = 4*M_PI;
    const double common_factor_p = (-four_pi/6)*(6*M_PI*(p01+p02)+(e01+e02))*pow(dr,2);
    const double common_factor_m = (four_pi/3)*pow(dr,3);

    // Start from the second order taylor series solution
    /*
    const double p1_2nd_order = p01 + (e01+p01)*common_factor_p;
    const double p2_2nd_order = p02 + (e02+p02)*common_factor_p;
    const double m1_2nd_order = e01*common_factor_m;
    const double m2_2nd_order = e02*common_factor_m;
    */
    // Should also calculate 2nd order solution to y_r,
    // but for now I'll just try without it and see what gives
    // Initialize state vector
    // tov_state y;
    /*
    y.push_back(m1_2nd_order);
    y.push_back(m2_2nd_order);
    y.push_back(p1_2nd_order);
    y.push_back(p2_2nd_order);
    y.push_back(y_0r);
    */
    tov_state y = {
        e01 * common_factor_m,
        e02 * common_factor_m,
        p01 + (e01 + p01) * common_factor_p,
        p02 + (e02 + p02) * common_factor_p,
        y_0r
    };
    /*
    y[0] = m1_2nd_order;
    y[1] = m2_2nd_order;
    y[2] = p1_2nd_order;
    y[3] = p2_2nd_order;
    y[4] = y_0r;
    */
    double r = dr;
    // Integrator
    const int sys_dim = 5;
    gsl_odeiv2_system sys = {twofluid_tov_eqns_gsl, NULL, sys_dim, reinterpret_cast<void *>(::std::addressof(*this))};
    const gsl_odeiv2_step_type* stepper_type = gsl_odeiv2_step_rkf45;
    gsl_odeiv2_step* stepper = gsl_odeiv2_step_alloc(stepper_type, sys_dim);
    gsl_odeiv2_control* stepper_control = gsl_odeiv2_control_y_new(0.0, 1e-10); // 1e-8 is the absolute error
    gsl_odeiv2_evolve* ode_ev = gsl_odeiv2_evolve_alloc(sys_dim);
    //controlled_runge_kutta<runge_kutta_dopri5<tov_state>> c_rk;
    //rosenbrock4_controller<double> c_rk;
    //controlled_runge_kutta<runge_kutta_fehlberg78<tov_state>> c_rk;
    // Control variables
    int i = 0;
    const int max_steps = 10000;
    const double r_max = 1e6;
    while ((y[2] > eos1.p_surface ||  y[3] > eos2.p_surface) && (i < max_steps)) {
        int status = gsl_odeiv2_evolve_apply(ode_ev, stepper_control, stepper, &sys, &r, r_max, &dr, &y[0]);
        // Add checks for NaN or negative/absurd values
        if (status != GSL_SUCCESS ||
            std::isnan(y[0]) || std::isnan(y[1]) || std::isnan(y[2]) || std::isnan(y[3]) || std::isnan(y[4]) ||
            y[0] < 0 || y[1] < 0 || r < 0 || y[2] < 0 || y[3] < 0) {
            std::cerr << "[WARN] Integration failed or unphysical values encountered at step " << i << std::endl;
            std::cerr << "[INFO] r: " << r << ", y[0]: " << y[0] << ", y[1]: " << y[1];
            std::cerr << ", y[2]: " << y[2] << ", y[3]: " << y[3] << ", y[4]: " << y[4] << std::endl;
            // Reset state to avoid further issues
            reset_state();
            break;
        }
        if (y[2] == 0.0 && y[3] == 0.0) {break;}
        i++;
    };
    if (y[2] <= eos1.p_surface && y[2] != 0) {y[2] = 0;};
    if (y[3] <= eos2.p_surface && y[3] != 0) {y[3] = 0;};
    y_r = y[4];
    if (i == max_steps) {
        std::cout << "Max number of steps reached! Solution may be incomplete ";
        std::cout << "[INFO] p1: " << y[2];
        std::cout << " p2: " << y[3]  << std::endl;
    };
    M = M_B + M_D;
    R = std::max(R_D,R_B);
    if (M_B != 0 && M_D != 0) {F_chi = M_D / M;}
    else {
        if (M_D == 0 && e02 != 0) {
            std::cout << "Surface of Baryon or DM fluid not reached" << std::endl;
            std::cout << "[INFO] e1: " << e01*CONVERSION::geom_to_dens_GeV_fm3;
            std::cout << " e2: " << e02*CONVERSION::geom_to_dens_GeV_fm3;
        }
        else {//e1=1.73799 e2=0.4137
            F_chi = 0.0;
            R_D = 0.0;
        }
    }
    calc_k2();
    calc_lambda();
    TOV_result result = {.M_B=M_B, .M_D=M_D,
                         .R_B=R_B, .R_D=R_D,
                         .M=M, .R=R,
                         .F_chi=F_chi, .k2=k2,
                         .lambda=lambda};
    gsl_odeiv2_evolve_free(ode_ev);
    gsl_odeiv2_control_free(stepper_control);
    gsl_odeiv2_step_free(stepper);
    return result;
};

void TwoFluid_TOV::calc_k2() {
    if (M > 0 && R > 0) {
        const double C = M / R;
        // Avoid unphysical compactness
        if (C >= 0.5) {
            k2 = 0;
            std::cerr << "[WARN] Compactness C >= 0.5 (C=" << C << "), setting k2=0" << std::endl;
            return;
        }
        const double num = ((8.0 / 5.0) * pow(1.0 - 2.0 * C, 2.0) *
                            pow(C, 5.0) * (2.0 * C * (y_r - 1.0) - y_r + 2.0));
        const double denom = (2.0 * C * (4.0 * (y_r + 1.0) * pow(C, 4.0) + (6.0 * y_r - 4.0) * pow(C, 3.0) +
                 (26.0 - 22.0 * y_r) * pow(C, 2.0) + 3.0 * (5.0 * y_r - 8.0) *
                 C - 3.0 * y_r + 6.0) - 3.0 * pow(1.0 - 2.0 * C, 2.0) *
                 (2.0 * C * (y_r - 1.0) - y_r + 2.0) *
                 log(1.0 / (1.0 - 2.0 * C)));
        if (denom <= 0) {
            k2 = 0;
            std::cerr << "[WARN] k2 denominator <= 0 (denom=" << denom << "), setting k2=0" << std::endl;
        } else {
            k2 = num / denom;
        }
    } else {
        k2 = 0;
        if (M <= 0 || R <= 0)
            std::cerr << "[WARN] Non-positive mass or radius in calc_k2 (M=" << M << ", R=" << R << "), setting k2=0" << std::endl;
    }
}

void TwoFluid_TOV::calc_lambda() {
    if (k2 > 0 && M > 0 && R > 0) {
        const double C = M / R;
        lambda = (2.0 / 3.0) * k2 / pow(C, 5);
        if (lambda < 0) {
            std::cerr << "[WARN] Negative lambda computed (" << lambda << "), setting lambda=0" << std::endl;
            lambda = 0;
        }
    } else {
        lambda = 0;
        if (k2 <= 0)
            std::cerr << "[WARN] k2 <= 0 in calc_lambda, setting lambda=0" << std::endl;
    }
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
    std::cout << "F_chi = " << F_chi << std::endl;
};

