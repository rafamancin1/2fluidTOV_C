#include "../include/eos_analytic.hpp"
#include "../include/eos_SIDM.hpp"
#include "../include/conversions.hpp"

EOS_SIDM::EOS_SIDM(const double m_chi, const double lambda)
    : EOS_Analytic("Self-Interacting Bosonic Dark Matter"), m_chi(m_chi), lambda(lambda), p_surface(1e-20) {}

EOS_SIDM::EOS_SIDM(const EOS_SIDM& other)
    : EOS_SIDM(other.m_chi, other.lambda) {}

EOS_SIDM::EOS_SIDM()
    : EOS_SIDM(0.0, 0.0) {}

double EOS_SIDM::pc_from_ec(const double& central_energy) const {
    const double central_energy_MeV = central_energy * CONVERSION::geom_to_dens_MeV4;
    const double first_factor = std::pow(m_chi, 4.0) / (9.0 * lambda);
    const double sqrt_term = std::sqrt(1.0 + (3.0 * lambda * central_energy_MeV) / std::pow(m_chi, 4.0));
    const double parenthesis_term = std::pow(sqrt_term - 1.0, 2.0);
    const double pressure = first_factor * parenthesis_term;
    const double pressure_geom = pressure * CONVERSION::pres_MeV4_to_geom;
    return pressure_geom;
}

double EOS_SIDM::energy_from_pressure(const double& pressure) const {
    const double pressure_MeV = pressure * CONVERSION::geom_to_pres_MeV4;
    const double first_term = 3.0 * pressure_MeV;
    const double second_term = 2.0 * std::pow(m_chi, 2.0) * std::sqrt(pressure_MeV / lambda);
    const double energy = first_term + second_term;
    const double energy_geom = energy * CONVERSION::dens_MeV4_to_geom;
    return energy_geom;
}

double EOS_SIDM::dedp(const double& pressure) const {
    if (pressure <= 0) {
        return 0.0;
    }
    const double pressure_MeV = pressure * CONVERSION::geom_to_pres_MeV4;
    const double first_term = 3.0;
    const double second_term_first_factor = std::pow(m_chi, 2.0) / lambda;
    const double secont_term_second_factor = std::pow(pressure_MeV / lambda, -0.5);
    const double second_term = second_term_first_factor * secont_term_second_factor;
    const double dedp_MeV = first_term + second_term;
    const double dedp_geom = dedp_MeV * CONVERSION::dens_MeV4_to_geom / CONVERSION::pres_MeV4_to_geom;
    return dedp_geom;
}
