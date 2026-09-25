#include "../include/eos_tabular.hpp"
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <limits>
#include <stdexcept>
//#include "../boost/math/interpolators/pchip.hpp"
//#include <gsl/gsl_errno.h>
//#include <gsl/gsl_spline.h>
#include <gsl/gsl_vector.h>


//using boost::math::interpolators::pchip; //Piecewise Cubic Hermite interpolation

// Directory holding the LALSimNeutronStarEOS_*.dat tables: $TWOFLUID_EOS_DIR if set,
// else the path baked in at build time (Makefile / setup.py), else ./eos_tables
static std::filesystem::path eos_table_dir() {
    const char* env_dir = std::getenv("TWOFLUID_EOS_DIR");
    if (env_dir != nullptr && env_dir[0] != '\0') {
        return env_dir;
    }
#ifdef TWOFLUID_DEFAULT_EOS_DIR
    return TWOFLUID_DEFAULT_EOS_DIR;
#else
    return "eos_tables";
#endif
}

//Copy constructor
EOS_Tabular::EOS_Tabular(const EOS_Tabular& other) : EOS_Tabular{other.eos_name} {};

EOS_Tabular::EOS_Tabular(std::string eos_name) : eos_name(eos_name), p_surface(0.0), e_min(0.0), e_max(0.0), p_max(0.0), tab_size(0) {
    const std::filesystem::path full_filename = eos_table_dir() / (eos_prefix + eos_name + eos_suffix);
    std::cout << full_filename.string() << std::endl;
    std::ifstream f_eos;
    f_eos.open(full_filename);
    if (!f_eos.is_open()) {
        throw std::runtime_error("Cannot open EOS table " + full_filename.string() +
                                 " (set TWOFLUID_EOS_DIR to the eos_tables directory)");
    };
    f_eos.seekg(0);
    std::vector<double> log_e_tab;
    std::vector<double> log_p_tab;
    //std::vector<double> log_e_tab_copy;
    //std::vector<double> log_p_tab_copy;
    //std::cout << "Reading file" << std::endl;
    while(true) {
        f_eos >> p;
        f_eos >> e;
        if (f_eos.eof()) { break; }
        // Tables start with a "0 0" row, whose log is -inf; it would poison the first spline intervals
        if (p <= 0 || e <= 0) { continue; }
        log_e_tab.push_back(log10(e));
        log_p_tab.push_back(log10(p));
        //log_e_tab_copy.push_back(log10(e));
        //log_p_tab_copy.push_back(log10(p));
    }
    tab_size = log_e_tab.size();
    if (tab_size < 3) {
        throw std::runtime_error("EOS table " + full_filename.string() + " has fewer than 3 usable rows");
    }
    //p_surface = pow(10, log_p_tab[0]);
    e_min = pow(10, log_e_tab[0]);
    e_max = pow(10, log_e_tab.back());
    p_max = pow(10, log_p_tab.back());
    p_surface = 1e-10*p_max;
    //p_surface = 1e-20;
    f_eos.close();
    initialize_splines(&log_e_tab[0], &log_p_tab[0]);
    //double log_e_tab_c[tab_size];
    //double log_p_tab_c[tab_size];
    /*
    for (int i = 0; i < tab_size; i++) {
        log_e_tab_c[i] = log_e_tab[i];
        log_p_tab_c[i] = log_p_tab[i];
    }
    */
    //std::vector<double> log_e_tab_copy = log_e_tab;
    //std::vector<double> log_p_tab_copy = log_p_tab;
    //e_p = pchip(std::move(log_p_tab_copy), std::move(log_e_tab_copy));
    //p_e = pchip(std::move(log_e_tab), std::move(log_p_tab));
    //tk::spline e_of_p(log_p_tab, log_e_tab);
    //tk::spline p_of_e(log_e_tab, log_p_tab);
    //gsl_interp_accel *acc_e = gsl_interp_accel_alloc();
    //gsl_interp_accel *acc_p = gsl_interp_accel_alloc();
    //gsl_spline *e_of_p = gsl_spline_alloc(gsl_interp_cspline, tab_size);
    //gsl_spline_init(e_of_p, log_e_tab_c, log_p_tab_c, tab_size);
    //e_p = tk::spline(log_p_tab, log_e_tab); //e_of_p;
    //p_e = tk::spline(log_e_tab, log_p_tab); //p_of_e;
    //e_p.setData(log_p_tab, log_e_tab);
    //p_e.setData(log_e_tab, log_p_tab);
    //free(log_e_tab_c);
    //free(log_p_tab_c);
};

EOS_Tabular::~EOS_Tabular() {
    gsl_spline_free(e_of_p);
    gsl_spline_free(p_of_e);
    gsl_interp_accel_free(acc_e);
    gsl_interp_accel_free(acc_p);
    //free(log_e_tab_c);
    //free(log_p_tab_c);
    //log_e_tab_c = nullptr;
    //log_p_tab_c = nullptr;
}

void EOS_Tabular::initialize_splines(double* log_e_tab_ptr, double* log_p_tab_ptr) {
    acc_e = gsl_interp_accel_alloc();
    acc_p = gsl_interp_accel_alloc();
    e_of_p = gsl_spline_alloc(gsl_interp_steffen, tab_size);
    p_of_e = gsl_spline_alloc(gsl_interp_steffen, tab_size);
    gsl_spline_init(e_of_p, log_p_tab_ptr, log_e_tab_ptr, tab_size);
    gsl_spline_init(p_of_e, log_e_tab_ptr, log_p_tab_ptr, tab_size);
    //free(log_e_tab_c);
    //free(log_p_tab_c);
}

// True when x lies inside the spline's tabulated range (false for NaN). Checking before
// evaluating keeps GSL from raising an error, whose default handler aborts the process.
static bool in_table(const gsl_spline* spline, double x) {
    return x >= spline->x[0] && x <= spline->x[spline->size - 1];
}

// log10 energy from log10 pressure; NaN outside the table
double EOS_Tabular::e_p(double log_pressure) {
    if (!in_table(e_of_p, log_pressure)) { return std::numeric_limits<double>::quiet_NaN(); }
    return gsl_spline_eval(e_of_p, log_pressure, acc_e);
}

// log10 pressure from log10 energy; NaN outside the table
double EOS_Tabular::p_e(double log_energy) {
    if (!in_table(p_of_e, log_energy)) { return std::numeric_limits<double>::quiet_NaN(); }
    return gsl_spline_eval(p_of_e, log_energy, acc_p);
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
    if (p_tab < p_max) {
        p_surface = 1e-10*p_tab; 
    }
    return p_tab;
};

// Analytic derivative of the log-log spline: e = 10^f(log10 p)  =>  de/dp = (e/p) f'(log10 p).
// Zero below the surface pressure, where energy_from_pressure is zero too; NaN above the table.
double EOS_Tabular::dedp(const double& pressure) {
    if (pressure <= 0 || pressure < p_surface) {
        return 0;
    }
    const double log_pressure = log10(pressure);
    if (!in_table(e_of_p, log_pressure)) { return std::numeric_limits<double>::quiet_NaN(); }
    const double energy = pow(10.0, e_p(log_pressure));
    const double dloge_dlogp = gsl_spline_eval_deriv(e_of_p, log_pressure, acc_e);
    return energy / pressure * dloge_dlogp;
};