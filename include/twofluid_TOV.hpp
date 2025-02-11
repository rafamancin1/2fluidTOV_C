#ifndef TWOFLUID_TOV_HPP
#define TWOFLUID_TOV_HPP
#include "eos.hpp"
#include "eos_tabular.hpp"
#include "eos_poly.hpp"
//#include "TOV_family.hpp"
#include "../boost/numeric/odeint.hpp"
#include <gsl/gsl_errno.h>

using namespace boost::numeric::odeint;

typedef std::vector<double> tov_state;
//typedef boost::numeric::ublas::vector<double> tov_state;
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
        // TwoFluid_TOV(EOS*, EOS*);
        // TwoFluid_TOV(EOS_Tabular*, EOS_Tabular*);
        TwoFluid_TOV(EOS_Tabular&, EOS_Poly&);
        // TwoFluid_TOV(EOS_Poly*, EOS_Poly*);    
        TOV_result integrate_two_fluid_tov(double, double);
        void print_result(TOV_result&);
        void reset_state();
        // public EOSs are a temporary solution
        EOS_Tabular& eos1;
        EOS_Poly& eos2;
        static int twofluid_tov_eqns_gsl(double const r, const double* const y, double* const dydr, void* const opaque) {
            assert(opaque);
            return(static_cast<TwoFluid_TOV *>(opaque)->twofluid_tov_eqns(r, y, dydr));
        }
    private:
        int twofluid_tov_eqns(double const, const double* const, double* const);
        //void tov_step(controlled_runge_kutta<runge_kutta_dopri5<tov_state>>&, tov_state&, double&, double&);
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
        const int y_0r = 2;
};

#endif