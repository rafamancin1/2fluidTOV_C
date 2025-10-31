import twofluidTOV
import argparse
from math import log10
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from tqdm import tqdm
#from bilby.gw.eos import EOSFamily, TabularEOS
N_SAMPLE = 200
mass_geom_to_Msun = 0.0006772199944005382
E1_0 = 2.0  # GeV fm^-3, this is the value of epsilon_B for the SLY EOS
EOS1_NAME = None
FCHI_PTS = []
MMAX_PTS = []
FCHI2_PTS = []
MMAX2_PTS = []
COLOR_PALETTE = plt.cm.rainbow(np.linspace(0, 1, 10))
COLOR_INDEX = 0

def main():
    parser = argparse.ArgumentParser(description="Setup plot parameters")
    parser.add_argument('--eos1', dest='eos1', nargs='+', required=True)
    parser.add_argument('--eos2', dest='eos2', required=False)
    parser.add_argument('--Achi', dest='Achi', nargs='+', type=float, required=False)
    parser.add_argument('--Gamma', dest='Gamma', type=float, required=False)
    parser.add_argument('--m_chi', dest='m_chi', nargs='+', type=float, required=False, default=0.1)
    parser.add_argument('--lambda_chi', dest='lambda_chi', nargs='+', type=float, required=False, default=0.1)
    parser.add_argument('--Fchi', dest='Fchi', nargs='+', type=float, required=False)
    parser.add_argument('--mode', dest='mode', required=False)
    args = parser.parse_args()
    plt.rcParams.update({'font.size': 24})
    fig, ax = plt.subplots(1,1,figsize=(12,8))
    if len(args.eos1) == 1:
        plot_name = args.eos1[0]+'_'
        global EOS1_NAME
        EOS1_NAME = args.eos1[0]
    else:
        plot_name = "multi_eos_"
    if args.mode == 'MR':
        plot_name += "MR_plots"
        ax.set_xlabel("R (km)")
        ax.set_ylabel(r"M ($M_\odot$)")
        #ax.set_xlim(8, 20)
    elif args.mode == 'lambdaM' or args.mode == 'lambdaM_test':
        plot_name += args.mode + "_plots"
        ax.set_xlabel(r'M ($M_\odot$)')
        ax.set_ylabel(r'$\Lambda$')
        ax.set_ylim(1e1, 1e4)
        #ax.set_xlim(1.0, 2.6)
        ax.set_yscale("log")
    elif args.mode == 'lambdae1':
        plot_name += "lambdae1_plots"
        ax.set_xlabel(r'$\epsilon_1$ (GeV fm$^{{-3}}$)')
        ax.set_ylabel(r'$\Lambda$')
        ax.set_xlim(0.53, 2.4)
        ax.set_ylim(1e1, 1e3)
    elif args.mode == 'MmaxFchi':
        plot_name += "MmaxFchi_plots"
        ax.set_xlabel(r'$F_\chi$')
        ax.set_ylabel(r'M$_{max}$ ($M_\odot$)')
    elif args.mode == 'lambdaFchi':
        plot_name += "lambdaFchi_plots"
        ax.set_xlabel(r'$F_\chi$')
        ax.set_ylabel(r'$\Lambda$')
    elif args.mode == 'e1M':
        plot_name += "e1_plots"
        ax.set_xlabel(r'$\epsilon_B$ (GeV fm$^{{-3}}$)')
        ax.set_ylabel(r'M ($M_\odot$)')
        #ax.set_xlim(0.53, 0.6)
    elif args.mode == 'e1R':
        plot_name += "e1R_plots"
        ax.set_xlabel(r'$\epsilon_B$ (GeV fm$^{{-3}}$)')
        ax.set_ylabel("R (km)")
    elif args.mode == 'e2M':
        plot_name += "e2_plots"
        ax.set_xlabel(r'$\epsilon_{DM}$ (GeV fm$^{{-3}}$)')
        ax.set_ylabel(r'M ($M_\odot$)')
    elif args.mode == 'e2R':
        plot_name += "e2R_plots"
        ax.set_xlabel(r'$\epsilon_{DM}$ (GeV fm$^{{-3}}$)')
        ax.set_ylabel("R (km)")
    elif args.mode == 'RM':
        plot_name += "RM_plots"
        ax.set_xlabel(r'M ($M_\odot$)')
        ax.set_ylabel("R (km)")
        #ax.set_xlim(1.0, 2.6)
        #ax.set_ylim
    elif args.mode == 'RFchi':
        plot_name += "RFchi_plots"
        ax.set_xlabel(r'$F_\chi$')
        ax.set_ylabel("R (km)")
    elif args.mode == "e1e2":
        plot_name += "e1e2_plots"
        ax.set_xlabel(r'$\epsilon_B$ (GeV fm$^{{-3}}$)')
        ax.set_ylabel(r'$\epsilon_{DM}$ (GeV fm$^{{-3}}$)')
    elif args.mode == 'test':
        eos1 = twofluidTOV.EOS_Tabular(args.eos1[0])
        eos2 = twofluidTOV.EOS_Poly(args.eos2, args.Achi[0], args.Gamma)
        #e1 = np.linspace(0.15, 2.4, 200) #BEST for SLY EOS
        print("Creating family")
        fam = twofluidTOV.TOV_Family(eos1, eos2)
        '''
        for e01 in e1:
            res = fam.calc_lambda_and_mass_directly(e01, args.Fchi[0], args.Achi[0])
            print("Mass: {}, Lambda: {}".format(res.M, res.lambda_param))
        '''
        #res = fam.calc_lambda_and_mass_directly(1.77, 0.34, 1e8)
        print(f"Testing ec1 prior with F_chi={args.Fchi[0]}")
        e1 = np.linspace(0.47, 1.0, 1000)
        res = [fam.calc_lambda_and_mass_directly(args.Fchi[0]) for e in e1]
        #print("Mass: {}, Lambda: {}".format(res.M * mass_geom_to_Msun, res.lambda_param))
        exit()
            
    import time
    start_time = time.time()
    if args.eos2 != None:
        eos1 = twofluidTOV.EOS_Tabular(args.eos1[0])
        if args.eos2 == 'Polytropic' and len(args.Achi) >= 1 and args.Gamma != None:
            eos2 = twofluidTOV.EOS_Poly(args.eos2, args.Achi[0], args.Gamma)
        elif args.eos2 == 'SIDM' and len(args.m_chi) >= 1 and len(args.lambda_chi) >= 1:
            eos2 = twofluidTOV.EOS_SIDM(args.m_chi[0], args.lambda_chi[0])
        else:
            print("Error: Invalid EOS2 type or parameters.")
            exit(1)
        if len(args.Fchi) > 1:
            for F_chi in args.Fchi:
                plot_label = r'$F_\chi = $' + str(int(F_chi*100)) + "%"
                create_plot(eos1, plot_label, ax, eos2, F_chi, args.mode)
        elif args.Achi is not None:
            if len(args.Achi) > 1:
                for A_chi in args.Achi:
                    eos2 = twofluidTOV.EOS_Poly(args.eos2, A_chi, args.Gamma)
                    if A_chi != 0:
                        n, e = get_number_and_exponent(A_chi)
                        A_chi_str = r'${} \cdot 10^{{ {} }}$'.format(n,e)
                    else:
                        A_chi_str = '0'
                    plot_label = r'$A_\chi = $' + A_chi_str
                    create_plot(eos1, plot_label, ax, eos2, args.Fchi[0], args.mode,A_chi) 
            else:
                create_plot(eos1, None, ax, eos2, args.Fchi[0], args.mode)
        elif args.m_chi is not None:
            if len(args.m_chi) > 1:
                for m_chi in args.m_chi:
                    eos2 = twofluidTOV.EOS_SIDM(m_chi, args.lambda_chi[0])
                    plot_label = r'$m_\chi = $' + str(m_chi) + " MeV"
                    create_plot(eos1, plot_label, ax, eos2, args.Fchi[0], args.mode, m_chi=m_chi)
                global FCHI_PTS
                global MMAX_PTS
                global FCHI2_PTS
                global MMAX2_PTS
                ax.scatter(FCHI_PTS, MMAX_PTS, color='red', s=100, label='Turning points')
                ax.scatter(FCHI2_PTS, MMAX2_PTS, color='red', s=100)
                ax.axhspan(2.5, 2.7, color='gray', hatch=r'/', alpha=0.3)
                ax.text(0.2, 2.75, 'GW190814', fontsize=20, color='black')
            else:
                create_plot(eos1, None, ax, eos2, args.Fchi[0], args.mode)
        elif args.lambda_chi is not None:
            if len(args.lambda_chi) > 1:
                for lambda_chi in args.lambda_chi:
                    eos2 = twofluidTOV.EOS_SIDM(args.m_chi[0], lambda_chi)
                    plot_label = r'$\lambda_\chi = $' + str(lambda_chi)
                    create_plot(eos1, plot_label, ax, eos2, args.Fchi[0], args.mode)
            else:
                create_plot(eos1, None, ax, eos2, args.Fchi[0], args.mode)
        else:
            #eos2 = twofluidTOV.EOS_Poly(args.eos2, args.Achi[0], args.Gamma)
            create_plot(eos1, None, ax, eos2, args.Fchi[0], args.mode)
    else:
        for eos_name in args.eos1:
            eos1 = twofluidTOV.EOS_Tabular(eos_name)
            create_plot(eos1, eos_name,  ax, mode=args.mode)
    end_time = time.time()
    print("Time taken: {} seconds".format(end_time - start_time))
    plot_name += ".pdf"
    #ax.set_title()
    fig.legend(loc=7)
    fig.tight_layout()
    fig.subplots_adjust(right=0.65)
    fig.savefig(plot_name)
    fig.show()

def create_plot(eos1, plot_label, ax, eos2=None, F_chi=None, mode='MR', A_chi=None, m_chi=None):
    global COLOR_INDEX
    global COLOR_PALETTE
    if eos2 != None:
        fam = twofluidTOV.TOV_Family(eos1, eos2, F_chi)
    else:
        fam = twofluidTOV.TOV_Family(eos1)
    if mode == 'MR':
        Ms = np.array(fam.Ms) * mass_geom_to_Msun
        Rs = np.array(fam.Rs) / 1e3
        #ax.set_xscale("log")
        ax.plot(Rs, Ms, label=plot_label)
    elif mode == 'lambdaM':
        Ms = np.array(fam.Ms) * mass_geom_to_Msun
        Ls = np.array(fam.lambdas)
        ax.plot(Ms, Ls, label=plot_label)
    elif mode == 'lambdaM_test':
        Ms = np.linspace(1.0, 2.6, 400)
        lambdas = np.array([fam.lambda_from_mass(M) for M in Ms])
        ax.plot(Ms, lambdas, label="lambda from spline")
        #ax.plot(np.array(fam.Ms)*mass_geom_to_Msun, np.array(fam.lambdas), label="lambda from family")
    elif mode == 'lambdae1':
        e1s = np.array(fam.e1s)
        Ls = np.array(fam.lambdas)
        ax.plot(e1s, Ls, label=plot_label)
    elif mode == 'MmaxFchi':
        Fchis = np.linspace(0.0, 0.95, 200)
        Mmaxs = np.zeros_like(Fchis)
        global EOS1_NAME
        global FCHI_PTS
        global MMAX_PTS
        global FCHI2_PTS
        global MMAX2_PTS
        data_fname = f'{EOS1_NAME}_{int(m_chi)}MeV_Mmax_Fchi.csv'
        data_fname2 = f'APR4_EPP_{int(m_chi)}MeV_Mmax_Fchi.csv'
        try:
            df = pd.read_csv(data_fname)
            print(f"File {data_fname} found, loading data.")
            Fchis = df['F_chi'].values
            Mmaxs = df['Mmax'].values
            dMdF = np.gradient(Mmaxs, Fchis)
            didx = np.argmin(np.abs(dMdF))
            Mmaxs_threshold = Mmaxs[didx]
            Fchis_threshold = Fchis[didx]

            # TEMPORARY CODE FOR APR4 EPP #
            df2 = pd.read_csv(data_fname2)
            Fchis2 = df2['F_chi'].values
            Mmaxs2 = df2['Mmax'].values
            dMdF2 = np.gradient(Mmaxs2, Fchis2)
            didx2 = np.argmin(np.abs(dMdF2))
            Mmaxs_threshold2 = Mmaxs2[didx2]
            Fchis_threshold2 = Fchis2[didx2]
            FCHI2_PTS.append(Fchis_threshold2)
            MMAX2_PTS.append(Mmaxs_threshold2)
            ###############################

            FCHI_PTS.append(Fchis_threshold)
            MMAX_PTS.append(Mmaxs_threshold)

            print(f"With m_chi= {m_chi} MeV, Threshold Mmax: {Mmaxs_threshold} at F_chi: {Fchis_threshold}")
            ax.plot(Fchis, Mmaxs, color=COLOR_PALETTE[COLOR_INDEX], label=plot_label)
            ax.plot(Fchis2, Mmaxs2, color=COLOR_PALETTE[COLOR_INDEX], linestyle='--')
            #ax.scatter(Fchis_threshold, Mmaxs_threshold, color='red', s=100, label='Transition point')
        except FileNotFoundError:
            print(f"File {data_fname} not found.")
            for i, f in tqdm(enumerate(Fchis), desc="Generating Mmax points"):
                fam = twofluidTOV.TOV_Family(eos1, eos2, f)
                Mmaxs[i] = max(fam.Ms) * mass_geom_to_Msun
            ax.plot(Fchis, Mmaxs, label=plot_label)
        df = pd.DataFrame({'F_chi': Fchis, 'Mmax': Mmaxs})
        df.to_csv(data_fname, index=False)
    elif mode == 'lambdaFchi':
        e1_0 = E1_0
        fam = twofluidTOV.TOV_Family(eos1, eos2)
        Fchis = np.linspace(0.0, 0.95, 200)
        Ls = np.array([fam.calc_lambda_and_mass_directly(e1_0, f, A_chi).lambda_param for f in Fchis])
        ax.plot(Fchis, Ls, label=plot_label)
    elif mode == 'e1M':
        e1s = np.array(fam.e1s)
        Ms = np.array(fam.Ms) * mass_geom_to_Msun
        ax.plot(e1s, Ms, label=plot_label)
    elif mode == 'e1R':
        e1s = np.array(fam.e1s)
        Rs = np.array(fam.Rs) / 1e3
        ax.plot(e1s, Rs, label=plot_label)
    elif mode == 'e2M':
        e2s = np.array(fam.e2s)
        Ms = np.array(fam.Ms) * mass_geom_to_Msun
        ax.plot(e2s, Ms, label=plot_label)
    elif mode == 'e2R':
        e2s = np.array(fam.e2s)
        Rs = np.array(fam.Rs) / 1e3
        ax.plot(e2s, Rs, label=plot_label)
    elif mode == 'RFchi':
        e1_0 = E1_0
        fam = twofluidTOV.TOV_Family(eos1, eos2)
        Fchis = np.linspace(0.0, 0.95, 200)
        RBs = np.array([fam.calc_lambda_and_mass_directly(e1_0, f, A_chi).R_B / 1e3 for f in Fchis])
        RDs = np.array([fam.calc_lambda_and_mass_directly(e1_0, f, A_chi).R_D / 1e3 for f in Fchis])
        colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
        color_rb = np.random.choice(colors)
        ax.plot(Fchis, RBs, label=plot_label, color=color_rb)
        ax.plot(Fchis, RDs, color=color_rb, linestyle='--')
    elif mode == 'RM':
        Ms = np.array(fam.Ms) * mass_geom_to_Msun
        RBs = np.array(fam.RBs) / 1e3
        RDs = np.array(fam.RDs) / 1e3
        colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
        color_rb = np.random.choice(colors)
        ax.plot(Ms, RBs, label=plot_label, color=color_rb)
        ax.plot(Ms, RDs, color=color_rb, linestyle='--')
    elif mode == 'e1e2':
        e1s = np.array(fam.e1s)
        e2s = np.array(fam.e2s)
        ax.plot(e1s, e2s, label=plot_label)
    COLOR_INDEX += 1

def get_number_and_exponent(value):
    """
    Returns the number and exponent for a given value.
    """
    if value == 0:
        return 0, 0
    exponent = int(log10(value))
    number = value / (10 ** exponent)
    return number, exponent


if __name__ == "__main__":
    main()
