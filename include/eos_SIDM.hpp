#ifndef EOS_SIDM_HPP
#define EOS_SIDM_HPP

#include <string>
#include <cmath>
#include "eos_analytic.hpp"

// Parameters live in EOS_Analytic: m_chi (DM particle mass in MeV),
// lambda_chi (self-interaction coupling constant) and p_surface.
class EOS_SIDM : public EOS_Analytic {
    public:
        EOS_SIDM(const double, const double);
        EOS_SIDM(const EOS_SIDM&);
        EOS_SIDM();

    public:
        double pc_from_ec(const double&) const override; // central pressure from central energy
        double energy_from_pressure(const double&) const override;
        double dedp(const double&) const override;
};

#endif // EOS_SIDM_HPP
