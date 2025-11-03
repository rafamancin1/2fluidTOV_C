#ifndef EOS_SIDM_HPP
#define EOS_SIDM_HPP

#include <string>
#include <cmath>
#include "eos_analytic.hpp"

class EOS_SIDM : public EOS_Analytic {
    public:
        EOS_SIDM(const double, const double);
        EOS_SIDM(const EOS_SIDM&);
        EOS_SIDM();
        
    public:
        double m_chi; // mass of the dark matter particle in MeV
        double lambda; // self-interaction coupling constant
        double p_surface;
    
    public:
        double pc_from_ec(const double&) const override; // central pressure from central energy
        double energy_from_pressure(const double&) const override;
        double dedp(const double&) const override;
};

#endif // EOS_SIDM_HPP