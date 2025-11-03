#ifndef EOS_ANALYTIC_HPP
#define EOS_ANALYTIC_HPP

#include <string>

class EOS_Analytic {
public:
    EOS_Analytic(std::string eos_name);
    std::string eos_name;
    double p_surface;
    double K;
    double Gamma;
    double m_chi;
    double lambda_chi;
    virtual double energy_from_pressure(const double&) const = 0;
    virtual double pc_from_ec(const double&) const = 0;
    virtual double dedp(const double&) const = 0;
    virtual ~EOS_Analytic() = default;
};

#endif