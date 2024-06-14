#include "../include/conversions.hpp"
 // Utils

 const double CONVERSION::MeV_fm3_to_GeV_fm3 = 1e-3;
 const double CONVERSION::GeV_fm3_to_MeV_fm3 = 1e3;

 const double CONVERSION::mass_geom_to_Msun = 0.0006772199944005382;
 const double CONVERSION::Msun_to_mass_geom = 1 / mass_geom_to_Msun;

 const double CONVERSION::m_to_km = 1e-3;
 const double CONVERSION::km_to_m = 1 / m_to_km;

 // Pressure conversions
const double CONVERSION::MeV_fm3_to_pa = 1.6022e32;
const double CONVERSION::pa_to_MeV_fm3 = 1 / MeV_fm3_to_pa;

const double CONVERSION::geom_to_pa = 1.2102e44;
const double CONVERSION::pa_to_geom = 1 / geom_to_pa;

const double CONVERSION::pres_MeV_fm3_to_geom = MeV_fm3_to_pa * pa_to_geom;
const double CONVERSION::geom_to_pres_MeV_fm3 = 1 / pres_MeV_fm3_to_geom;

// Density conversions
const double CONVERSION::nucleon_m = 939.57; // in MeV
const double CONVERSION::nucleon_m_GeV = nucleon_m / 1e3;

const double CONVERSION::MeV_fm3_to_kg_m3 = 1.7827e15;
const double CONVERSION::kg_m3_to_MeV_fm3 = 1 / MeV_fm3_to_kg_m3;

const double CONVERSION::kg_m3_to_geom = 7.4261602691186655e-28;
const double CONVERSION::geom_to_kg_m3 = 1 / kg_m3_to_geom;

const double CONVERSION::dens_MeV_fm3_to_geom = MeV_fm3_to_kg_m3 * kg_m3_to_geom;
const double CONVERSION::geom_to_dens_MeV_fm3 = 1 / dens_MeV_fm3_to_geom;

const double CONVERSION::dens_GeV_fm3_to_geom = GeV_fm3_to_MeV_fm3 * dens_MeV_fm3_to_geom;
const double CONVERSION::geom_to_dens_GeV_fm3 = 1 / dens_GeV_fm3_to_geom;

const double CONVERSION::g_cm3_to_kg_m3 = 1e3;
const double CONVERSION::kg_m3_to_g_cm3 = 1 / g_cm3_to_kg_m3;

const double CONVERSION::fm3_to_MeV_fm3 = nucleon_m;
const double CONVERSION::MeV_fm3_to_fm3 = 1 / fm3_to_MeV_fm3;

const double CONVERSION::fm3_to_geom = fm3_to_MeV_fm3 * dens_MeV_fm3_to_geom;
const double CONVERSION::geom_to_fm3 = 1 / fm3_to_geom;
