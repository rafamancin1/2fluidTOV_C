#include "../include/eos_analytic.hpp"

EOS_Analytic::EOS_Analytic(std::string eos_name) : eos_name(eos_name) {
    // Default constructor for EOS_Analytic
    p_surface = 1e-20; // Default surface pressure
    K = 0.0; // Default value for K
    Gamma = 0.0; // Default value for Gamma
    m_chi = 0.0; // Default value for m_chi
    lambda_chi = 0.0; // Default value for lambda_chi
}


