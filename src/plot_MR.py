import twofluidTOV
import argparse
from math import log10
import numpy as np
import matplotlib.pyplot as plt
#from bilby.gw.eos import EOSFamily, TabularEOS
N_SAMPLE = 200
mass_geom_to_Msun = 0.0006772199944005382

def main():
    parser = argparse.ArgumentParser(description="Setup plot parameters")
    parser.add_argument('--eos1', dest='eos1', nargs='+', required=True)
    parser.add_argument('--eos2', dest='eos2', required=False)
    parser.add_argument('--Achi', dest='Achi', nargs='+', type=float, required=False)
    parser.add_argument('--Gamma', dest='Gamma', type=float, required=False)
    parser.add_argument('--Fchi', dest='Fchi', nargs='+', type=float, required=False)
    parser.add_argument('--mode', dest='mode', required=False)
    args = parser.parse_args()
    plt.rcParams.update({'font.size': 20})
    fig, ax = plt.subplots(1,1,figsize=(12,8))
    if args.mode == 'MR':
        plot_name = "MR_plots"
        ax.set_xlabel("R (km)")
        ax.set_ylabel(r"M ($M_\odot$)")
        ax.set_xlim(8, 20)
    elif args.mode == 'lambdaM':
        plot_name = "lambda_plots"
        ax.set_xlabel(r'M ($M_\odot$)')
        ax.set_ylabel(r'$\Lambda$')
        ax.set_ylim(1e1, 1e4)
        ax.set_xlim(1.0, 2.6)
        ax.set_yscale("log")
    elif args.mode == 'e1M':
        plot_name = "e1_plots"
        ax.set_xlabel(r'$\epsilon_1$')
        ax.set_ylabel(r'M ($M_\odot$)')
        ax.set_xlim(0.15, 2.4)

    elif args.mode == 'test':
        eos1 = twofluidTOV.EOS_Tabular(args.eos1[0])
        eos2 = twofluidTOV.EOS_Poly(args.eos2, args.Achi[0], args.Gamma)
        e1 = np.linspace(0.15, 2.4, 200) #BEST for SLY EOS
        fam = twofluidTOV.TOV_Family(eos1, eos2)
        for e01 in e1:
            res = fam.calc_lambda_and_mass_directly(e01, args.Fchi[0], args.Achi[0])
            print("Mass: {}, Lambda: {}".format(res.M, res.lambda_param))
        exit()
            
    
    if args.eos2 != None:
        eos1 = twofluidTOV.EOS_Tabular(args.eos1[0])
        if len(args.Fchi) > 1:
            for F_chi in args.Fchi:
                eos2 = twofluidTOV.EOS_Poly(args.eos2, args.Achi[0], args.Gamma)
                plot_label = r'$F_\chi = $' + str(int(F_chi*100)) + "%"
                create_plot(eos1, plot_label, ax, eos2, F_chi, args.mode)
        elif len(args.Achi) > 1:
            for A_chi in args.Achi:
                eos2 = twofluidTOV.EOS_Poly(args.eos2, A_chi, args.Gamma)
                if A_chi != 0:
                    A_chi_str = r'$10^{{ {} }}$'.format(int(log10(A_chi)))
                else:
                    A_chi_str = '0'
                plot_label = r'$A_\chi = $' + A_chi_str
                create_plot(eos1, plot_label, ax, eos2, args.Fchi[0], args.mode) 
    else:
        for eos_name in args.eos1:
            eos1 = twofluidTOV.EOS_Tabular(eos_name)
            create_plot(eos1, eos_name,  ax, mode=args.mode)

    plot_name += ".png"
    #ax.set_title()
    fig.legend()
    fig.savefig(plot_name)
    fig.show()

def create_plot(eos1, plot_label, ax, eos2=None, F_chi=None, mode='MR'):
    if eos2 != None:
        fam = twofluidTOV.TOV_Family(eos1, eos2, F_chi)
    else:
        fam = twofluidTOV.TOV_Family(eos1)
    if mode == 'MR':
        Ms = np.array(fam.Ms) * mass_geom_to_Msun
        Rs = np.array(fam.Rs) / 1e3
        ax.plot(Rs, Ms, label=plot_label)
    elif mode == 'lambdaM':
        Ms = np.array(fam.Ms) * mass_geom_to_Msun
        Ls = np.array(fam.lambdas)
        ax.plot(Ms, Ls, label=plot_label)
    elif mode == 'e1M':
        e1s = np.array(fam.e1s)
        Ms = np.array(fam.Ms) * mass_geom_to_Msun
        ax.plot(e1s, Ms, label=plot_label)

if __name__ == "__main__":
    main()