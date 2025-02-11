#include "../include/twofluid_TOV.hpp"
//#include "../boost/numeric/odeint.hpp"
#include "../include/conversions.hpp"
#include <boost/bind.hpp>
#include <cmath>
#include <gsl/gsl_odeiv2.h>
#include <numeric>




// TwoFluid_TOV::TwoFluid_TOV(EOS_Tabular* eos1, EOS_Tabular* eos2) : eos1(eos1), eos2(eos2){};
TwoFluid_TOV::TwoFluid_TOV(EOS_Tabular& eos1, EOS_Poly& eos2) : eos1(eos1), eos2(eos2), M_B(0), M_D(0), R_B(0), R_D(0), M(0), R(0), dr_min(1e-2) {};
// TwoFluid_TOV::TwoFluid_TOV(EOS_Poly* eos1, EOS_Poly* eos2) : eos1(eos1), eos2(eos2) {};

int TwoFluid_TOV::twofluid_tov_eqns(double const r, const double* const y, double* const dydr) {
    const double m1 = y[0];
    const double m2 = y[1];
    const double p1 = y[2] > eos1.p_surface ? y[2] : 0.0;
    const double p2 = y[3] > eos2.p_surface ? y[3] : 0.0;
    const double y_r = y[4];
    if (p1 == 0.0 && M_B == 0) {M_B = m1, R_B = r;}
    if (p2 == 0.0 && M_D == 0) {M_D = m2, R_D = r;}

    const double e1 = (eos1.energy_from_pressure)(p1);
    const double e2 = (eos2.energy_from_pressure)(p2);
    const double dedp1 = (eos1.dedp)(p1);
    const double dedp2 = (eos2.dedp)(p2);
    const double ep1 = dedp1*(e1 + p1);
    const double ep2 = dedp2*(e2 + p2);
    const double p_total = p1 + p2;
    const double e_total = e1 + e2;
    const double m = m1 + m2;
    const double sdenom = r - 2*m;
    const double denom = r*sdenom;
    const double num = (4*M_PI*(pow(r,3))*p_total + m);

    const double Q2nd = 4*pow(((m + 4*M_PI*pow(r,3)*p_total) / denom), 2);
    const double Q1st = 4*M_PI*r*(5*e_total + 9*p_total + ep1 + ep2 - 6/(4*M_PI*pow(r,2))) / sdenom;
    const double Qr = Q1st - Q2nd;
    const double Fr = (r - 4*M_PI*pow(r,3)*(e_total - p_total)) / sdenom;

    const double dm1dr = 4*M_PI*e1*pow(r,2);
    const double dm2dr = 4*M_PI*e2*pow(r,2);
    const double dp1dr = -(p1 + e1)*num / denom;
    const double dp2dr = -(p2 + e2)*num / denom;
    const double dy_rdr = -(pow(y_r,2) + y_r*Fr + pow(r,2)*Qr)/r;

    dydr[0] = dm1dr;
    dydr[1] = dm2dr;
    dydr[2] = dp1dr;
    dydr[3] = dp2dr;
    dydr[4] = dy_rdr;
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

    // Start from the second order taylor series solution
    const double p1_2nd_order = p01 + (-4*M_PI/6)*(e01+p01)*(6*M_PI*(p01+p02)+(e01+e02))*pow(dr,2);
    const double p2_2nd_order = p02 + (-4*M_PI/6)*(e02+p02)*(6*M_PI*(p01+p02)+(e01+e02))*pow(dr,2);
    const double m1_2nd_order = (4/3)*M_PI*e01*pow(dr,3);
    const double m2_2nd_order = (4/3)*M_PI*e02*pow(dr,3);
    // Should also calculate 2nd order solution to y_r,
    // but for now I'll just try without it and see what gives
    // Initialize state vector
    tov_state y;
    
    y.push_back(m1_2nd_order);
    y.push_back(m2_2nd_order);
    y.push_back(p1_2nd_order);
    y.push_back(p2_2nd_order);
    y.push_back(y_0r);
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
    const gsl_odeiv2_step_type* stepper_type = gsl_odeiv2_step_rk8pd;
    gsl_odeiv2_step* stepper = gsl_odeiv2_step_alloc(stepper_type, sys_dim);
    gsl_odeiv2_control* stepper_control = gsl_odeiv2_control_y_new(0.0, 1e-8);
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
        if (y[2] == 0.0 && y[3] == 0.0) {break;}
        i++;
    };
    if (y[2] <= eos1.p_surface && y[2] != 0) {y[2] = 0;};
    if (y[3] <= eos2.p_surface && y[3] != 0) {y[3] = 0;};
    double y_r = y[4];
    if (i == max_steps) {std::cout << "Max number of steps reached! Solution may be incomplete" << std::endl;};
    M = M_B + M_D;
    R = R_B > R_D ? R_B : R_D;
    if (M_B != 0 && M_D != 0) {F_chi = M_D / M;}
    else {
        if (M_D == 0 && e02 != 0) {
            std::cout << "Surface of Baryon or DM fluid not reached" << std::endl;
            std::cout << "[INFO] e1: " << e01*CONVERSION::geom_to_dens_GeV_fm3;
            std::cout << " e2: " << e02*CONVERSION::geom_to_dens_GeV_fm3;
            std::cout << " K: " << eos2.K;
            std::cout << " Gamma: " << eos2.Gamma << std::endl;  
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
    if (M != 0 && R != 0) {
        const double C = M/R;
        const double num = ((8.0 / 5.0) * pow(1.0 - 2.0 * C, 2.0) *
                            pow(C, 5.0) * (2.0 * C * (y_r - 1.0) - y_r + 2.0));
        const double denom = (2.0 * C * (4.0 * (y_r + 1.0) * pow(C, 4.0) + (6.0 * y_r - 4.0) * pow(C, 3.0) +
                 (26.0 - 22.0 * y_r) * pow(C, 2.0) + 3.0 * (5.0 * y_r - 8.0) *
                 C - 3.0 * y_r + 6.0) - 3.0 * pow(1.0 - 2.0 * C, 2.0) *
                 (2.0 * C * (y_r - 1.0) - y_r + 2.0) *
                 log(1.0 / (1.0 - 2.0 * C)));
        k2 = num / denom;
    }
    else {
        k2 = 0;
    }
};

void TwoFluid_TOV::calc_lambda() {
    if (k2 != 0) {
        const double C = M / R;
        lambda = (2.0 / 3.0) * k2 / pow(C, 5); 
    }
    else {
        lambda = 0;
    }
};

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

