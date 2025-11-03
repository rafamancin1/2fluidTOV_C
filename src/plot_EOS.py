import twofluidTOV
import numpy as np
import matplotlib.pyplot as plt
import argparse
from math import log10

N_SAMPLE = 200
mass_geom_to_Msun = 0.0006772199944005382
MeV_fm3_to_kg_m3 = 1.7827e15  # Conversion factor from MeV fm^-3 to kg m^-3
kg_m3_to_MeV_fm3 = 1 / MeV_fm3_to_kg_m3  # Conversion factor from kg m^-3 to MeV fm^-3
kg_m3_to_geom = 7.4261602691186655e-28  # Conversion factor from kg m^-3 to geometric units
geom_to_kg_m3 = 1 / kg_m3_to_geom  # Conversion factor from geometric units to kg m^-3
MeV_fm3_to_GeV_fm3 = 1 / 1000.0  # Conversion factor from MeV fm^-3 to GeV fm^-3
GeV_fm3_to_MeV_fm3 = 1000.0  # Conversion factor from GeV fm^-3 to MeV fm^-3
geom_to_pa = 1.2102e44  # Conversion factor from geometric units to pascals
pa_to_geom = 1 / geom_to_pa  # Conversion factor from pascals to geometric units
pa_to_dynes_cm2 = 10
kg_m3_to_g_cm3 = 1e-3  # Conversion factor from kg m^-3 to g cm^-3
MeV_fm3_to_kg_m3 = 1.7827e15  # Conversion factor from MeV fm^-3 to kg m^-3
kg_m3_to_MeV_fm3 = 1 / MeV_fm3_to_kg_m3  # Conversion factor from kg m^-3 to MeV fm^-3
dens_MeV_fm3_to_geom = MeV_fm3_to_kg_m3 * kg_m3_to_geom
dens_GeV_fm3_to_geom = GeV_fm3_to_MeV_fm3 * dens_MeV_fm3_to_geom
geom_to_dens_MeV_fm3 = 1 / dens_MeV_fm3_to_geom  # Conversion factor from geometric units to MeV fm^-3
geom_to_dens_GeV_fm3 = 1 / dens_GeV_fm3_to_geom  # Conversion factor from geometric units to GeV fm^-3
MeV_fm3_to_pa = 1.6022e32  # Conversion factor from MeV fm^-3 to pascals
pa_to_MeV_fm3 = 1 / MeV_fm3_to_pa  # Conversion factor from pascals to MeV fm^-3
pres_MeV_fm3_to_geom = MeV_fm3_to_pa * pa_to_geom
geom_to_pres_MeV_fm3 = 1 / pres_MeV_fm3_to_geom  # Conversion factor from geometric units to MeV fm^-3
#pres_MeV_fm3_to_geom = 1.6022e32 * 1/1.2102e44 # Conversion factor from MeV fm^-3 to geometric units
geom_to_g_cm3 = geom_to_kg_m3 * kg_m3_to_g_cm3  # Conversion factor from geometric units to g cm^-3
g_cm3_to_geom = 1 / geom_to_g_cm3  # Conversion factor from g cm^-3 to geometric units
geom_to_dynes_cm2 = geom_to_pa * pa_to_dynes_cm2  # Conversion factor from geometric units to dynes cm^-2
E1_0 = 2.0  # GeV fm^-3, this is the value of epsilon_B for the SLY EOS

def main():
    parser = argparse.ArgumentParser(description="Setup plot parameters")
    parser.add_argument('--eos', dest='eos_names', nargs='+', required=True)
    parser.add_argument('--Achi', dest='Achi', nargs='+', type=float, required=False)
    parser.add_argument('--Gamma', dest='Gamma', type=float, required=False)
    parser.add_argument('--m_chi', dest='m_chi', type=float, required=False, default=0.1)
    parser.add_argument('--lambda_chi', dest='lambda_chi', type=float, required=False, default=0.1)
    parser.add_argument('--plot_mode', dest='plot_mode', type=str, required=False, default='default',)
    parser.add_argument('--units', dest='units', type=str, required=False, default='cgs',)
    args = parser.parse_args()
    plt.rcParams.update({'font.size': 20})
    fig, ax = plt.subplots(1, 1, figsize=(12, 8))
    if len(args.eos_names) == 1:
        plot_name = args.eos_names[0] + '_'
    else:
        plot_name = "multi_eos_"
    plot_name += "EOS_plots"
    #ax.set_xlabel(r'$\epsilon_c$ (GeV fm$^{-3}$)')
    #ax.set_ylabel(r'$P_c$ (MeV fm$^{-3}$)')
    #ax.set_xlim(0.53, 2.0)
    #ax.set_ylim(0.0, 1.0)
    ax.set_yscale("log")
    ax.set_xscale("log")
    for eos_name in args.eos_names:
        print(f"Plotting EOS: {eos_name}")
        if eos_name == 'Polytropic':
            eos = twofluidTOV.EOS_Poly(eos_name, args.Achi[0], args.Gamma)
        elif eos_name == 'SIDM':
            eos = twofluidTOV.EOS_SIDM(args.m_chi, args.lambda_chi)
        else:
            eos = twofluidTOV.EOS_Tabular(eos_name)
        #More EOS types can be added here
        #pressure = np.linspace(0.01, 1.0, N_SAMPLE) * pres_MeV_fm3_to_geom  # Pressure in geometric units
        energy_density = np.linspace(0.03, 2.4, N_SAMPLE) * dens_GeV_fm3_to_geom  # Energy density in geometric units
        #print(energy_density)
        pressure = np.array([eos.pc_from_ec(energy) for energy in energy_density])
        if args.plot_mode == 'default':
            #energy_density = np.array([eos.energy_from_pressure(p) for p in pressure])
            #pressure = [eos.pc_from_ec(energy) for energy in energy_density]
            if args.units == 'cgs':
                ax.set_xlabel(r'$\epsilon_c$ (g cm$^{-3}$)')
                ax.set_ylabel(r'$P_c$ (dynes cm$^{-2}$)')
                energy_density *= geom_to_g_cm3  # Convert energy density to g cm^-3
                pressure *= geom_to_dynes_cm2  # Convert pressure to dynes cm^-2
            elif args.units == 'natural':
                ax.set_xlabel(r'$\epsilon_c$ (MeV fm$^{-3}$)')
                ax.set_ylabel(r'$P_c$ (MeV fm$^{-3}$)')
                energy_density *= geom_to_dens_MeV_fm3  # Convert energy density to GeV fm^-3
                pressure *= geom_to_pres_MeV_fm3
            ax.plot(energy_density, pressure, label=eos_name)
        elif args.plot_mode == 'dedp':
            # Plotting dP/dε vs ε
            dedp = np.array([eos.dedp(p) for p in pressure])
            if args.units == 'cgs':
                ax.set_xlabel(r'$P_c$ (dynes cm$^{-2}$)')
                ax.set_ylabel(r'$dP_c/d\epsilon_c$ (dynes cm$^{-2}$ g$^{-1}$ cm$^{3}$)')
                energy_density *= geom_to_g_cm3
                dedp *= geom_to_g_cm3 / geom_to_dynes_cm2  # Convert dP/dε to dynes cm^-2 g^-1 cm^3
            elif args.units == 'natural':
                ax.set_xlabel(r'$\epsilon_c$ (MeV fm$^{-3}$)')
                ax.set_ylabel(r'$dP_c/d\epsilon_c$ (MeV fm$^{-3}$ g$^{-1}$ cm$^{3}$)')
                energy_density *= geom_to_dens_MeV_fm3
                dedp *= geom_to_dens_MeV_fm3 / geom_to_pres_MeV_fm3
            ax.plot(energy_density, dedp, label=eos_name)
            plot_name += '_dedp'
        ax.legend()
    plt.savefig(plot_name + '.png')

if __name__ == "__main__":
    main()
