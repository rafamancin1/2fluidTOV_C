#include <iostream>
#include <stdio.h>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <ctime>
#include <exception>
#include <string>

#include "../include/eos_tabular.hpp"
#include "../include/eos_poly.hpp"
#include "../include/eos_SIDM.hpp"
#include "../include/twofluid_TOV.hpp"
#include "../include/TOV_family.hpp"
#include "../include/conversions.hpp"

std::string tabular_eos_name;
double e01 = 0.3;
double e02 = 0.8;
double m_chi = 1e7;
double lambda_chi = 2;
double F_chi = 0.2;
double e2_max = E2_MAX_DEFAULT;

static void print_usage(const char* prog) {
    std::cout << "Usage: " << prog << " -t <tabular_eos_name> -k <m_chi (MeV)> -g <lambda_chi>"
              << " -b <e01 (GeV/fm^3)> -d <e02 (GeV/fm^3)> -c <F_chi> [-e <e2_max (GeV/fm^3)>] [-f] [-T]" << std::endl;
}

static int run(int argc, char* argv[]) {
    bool family = false;
    bool test = false;
    for(int i=0; i<argc; i++) {
        if (argv[i][0] == '-') {
            const char opt = argv[i][1];
            const bool needs_value = (opt == 't' || opt == 'g' || opt == 'k' || opt == 'b' || opt == 'd' || opt == 'c' || opt == 'e');
            if (needs_value && i+1 >= argc) {
                std::cout << "Option " << argv[i] << " needs a value" << std::endl;
                print_usage(argv[0]);
                return 1;
            }
            switch(opt) {
                case 't':
                    tabular_eos_name = argv[i+1];
                    break;
                case 'g':
                    sscanf(argv[i+1], "%lf", &lambda_chi);
                    break;
                case 'k':
                    sscanf(argv[i+1], "%lf", &m_chi);
                    break;
                case 'f':
                    family = true;
                    break;
                case 'b':
                    sscanf(argv[i+1], "%lf", &e01);
                    break;
                case 'd':
                    sscanf(argv[i+1], "%lf", &e02);
                    break;
                case 'T':
                    test = true;
                    break;
                case 'c':
                    sscanf(argv[i+1], "%lf", &F_chi);
                    break;
                case 'e':
                    sscanf(argv[i+1], "%lf", &e2_max);
                    break;
                default:
                    std::cout << "Unknown option: " << argv[i] << std::endl;
                    print_usage(argv[0]);
                    return 1;
            }
        }
    }
    std::cout << "Starting" << std::endl;
    std::cout << "Tabular EOS name: " << tabular_eos_name << std::endl;
    EOS_Tabular eos1(tabular_eos_name);
    std::cout << "SIDM EOS has m_chi=" << m_chi << " and lambda_chi=" << lambda_chi << std::endl;
    EOS_SIDM eos2(m_chi, lambda_chi);
    TwoFluid_TOV model(eos1, eos2);
    if (test) {
        std::cout << "Testing" << std::endl;
        TOV_Family sols(eos1, eos2);
        sols.e2_max = e2_max;
        std::chrono::time_point time1 = std::chrono::high_resolution_clock::now();
        double lambda = sols.calc_lambda(F_chi, 1.2);
        std::chrono::time_point time2 = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(time2-time1);
        std::cout << "Execution took " << time_span.count() << " seconds" << std::endl;
        std::cout << "Lambda is " << lambda << std::endl;
        return 0;
    }
    if (family) {
        std::cout << "Beginning integration of family of solutions with F_chi=" << F_chi << std::endl;
        double time1 = (double) clock() / CLOCKS_PER_SEC;
        TOV_Family sols(eos1, eos2, F_chi, e2_max);
        double time2 = (double) clock() / CLOCKS_PER_SEC;
        std::cout << "Execution took " << time2-time1 << " seconds" << std::endl;
        double M_max = 0.0;
        int n_failed = 0;
        for (double M_i : sols.Ms) {
            if (std::isfinite(M_i)) { M_max = std::max(M_max, M_i); } else { n_failed++; }
        }
        if (n_failed > 0) {
            std::cout << n_failed << " of " << sols.Ms.size() << " points failed or need e2 > e2_max=" << e2_max
                      << " GeV/fm^3 (raise it with -e)" << std::endl;
        }
        std::cout << "Maximum mass(M_sun) = " << M_max * CONVERSION::mass_geom_to_Msun << std::endl;
        std::cout << "Lambda(1.4 M_sun) = " << sols.lambda_from_mass(1.4) << std::endl;
    }
    else {
        TOV_result result = model.integrate_two_fluid_tov(e01, e02);
        model.print_result(result);
    }
    return 0;
}

int main(int argc, char* argv[]) {
    try {
        return run(argc, argv);
    }
    catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
}
