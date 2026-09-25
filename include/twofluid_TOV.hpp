#ifndef TWOFLUID_TOV_HPP
#define TWOFLUID_TOV_HPP
#include "eos_analytic.hpp"
#include "eos_tabular.hpp"
#include "eos_poly.hpp"
#include <cassert>
#include <vector>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv2.h>

typedef std::vector<double> tov_state;
typedef struct TOV_result {
            const double M_B;
            const double M_D;
            const double R_B;
            const double R_D;
            const double M;
            const double R;
            const double F_chi;
            const double k2;
            const double lambda;
        } TOV_result;



class TwoFluid_TOV {
    public:
        TwoFluid_TOV(EOS_Tabular&, EOS_Analytic&);
        // Returns nan_result() when the configuration cannot be solved
        TOV_result integrate_two_fluid_tov(double, double);
        static TOV_result nan_result(); // every field NaN: signals a failed solve
        void print_result(TOV_result&);
        void reset_state();
        // public EOSs are a temporary solution
        EOS_Tabular& eos1;
        EOS_Analytic& eos2;
        static int twofluid_tov_eqns_gsl(double const r, const double* const y, double* const dydr, void* const opaque) {
            assert(opaque);
            return(static_cast<TwoFluid_TOV *>(opaque)->twofluid_tov_eqns(r, y, dydr));
        }
    private:
        int twofluid_tov_eqns(double const, const double* const, double* const);
        double locate_surface(int, double, const tov_state&, double, gsl_odeiv2_step*, const gsl_odeiv2_system*);
        void calc_k2();
        void calc_lambda();
        double M_B;
        double R_B;
        double M_D;
        double R_D;
        double M;
        double R;
        double y_r;
        double k2;
        double lambda;
        double F_chi;
        double dr_min;
        // Surface pressures and which fluids are still integrated, set per integration
        double ps1;
        double ps2;
        bool active1;
        bool active2;
        const double y_0r = 2.0;
};

#endif
