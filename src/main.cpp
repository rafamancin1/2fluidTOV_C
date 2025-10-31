#include <iostream>
#include <stdio.h>

#include "../include/eos_tabular.hpp"
#include "../include/eos_poly.hpp"
#include "../include/eos_SIDM.hpp"
#include "../include/twofluid_TOV.hpp"
#include "../include/TOV_family.hpp"
//#include "../include/utils.hpp"
#include <time.h>

char tabular_eos_name[10];
//std::string tabular_eos_name;
double e01 = 0.3;
double e02 = 0.8;
double K = 1e7;
double Gamma = 2;
double F_chi = 0.2;
int print;

int main(int argc, char* argv[]) {
    bool family = false;
    bool test = false;
    bool generation = false;
    for(int i=0; i<argc; i++) {
        if (argv[i][0] == '-') {
            switch(argv[i][1]) {
                case 't':
                    sscanf(argv[i+1], "%s", tabular_eos_name);
                    break;
                case 'g':
                    sscanf(argv[i+1], "%lf", &Gamma);
                    break;
                case 'k':
                    sscanf(argv[i+1], "%lf", &K);
                    break;
                case 'f':
                    family = true;
                    break;
                case 'O':
                    generation = true;
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
                default:
                    std::cout << "Unknown option: " << argv[i] << std::endl;
                    std::cout << "Usage: " << argv[0] << " -t <tabular_eos_name> -g <Gamma> -k <K> -b <e01> -d <e02> -c <F_chi> [-f] [-O] [-T]" << std::endl;
                    exit(1);
                    break;
            }
        }   
    }
    std::cout << "Starting" << std::endl;
    std::cout << "Tabular EOS name: " << tabular_eos_name << std::endl;
    EOS_Tabular eos1(tabular_eos_name);
    //std::cout << "Polytropic EOS has K=" << K << " and Gamma=" << Gamma << std::endl;
    std::cout << "SIDM EOS has m_chi=" << K << " and lambda_chi=" << Gamma << std::endl;
    EOS_SIDM eos2(K, Gamma);
    //EOS_Poly eos2("Polytropic", K, Gamma);
    TwoFluid_TOV model(eos1, eos2);
    const int n_samples = 200;
    if (test) {
        std::cout << "Testing" << std::endl;
        TOV_Family sols(eos1, eos2);
        // double time1 = (double) clock() / CLOCKS_PER_SEC;
        std::chrono::time_point time1 = std::chrono::high_resolution_clock::now();
        //double lambda = sols.calc_lambda_parallel(K, F_chi, 1.2, 10);
        double lambda = sols.calc_lambda(F_chi, 1.2);
        //std::cout << "e1 is " << e01 << std::endl;
        //test_calc_e2(K, 0.1, e01);
        //double time2 = (double) clock() / CLOCKS_PER_SEC;
        std::chrono::time_point time2 = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(time2-time1);
        std::cout << "Execution took " << time_span.count() << " seconds" << std::endl;
        std::cout << "Lambda is " << lambda << std::endl;
        return 0;
    }
    if (generation) {
        TOV_Family sols(eos1, eos2);
        double time1 = (double) clock() / CLOCKS_PER_SEC;
        sols.generate_lambda_f_points();
        double time2 = (double) clock() / CLOCKS_PER_SEC;
        std::cout << "Execution took " << time2-time1 << " seconds" << std::endl;
        return 0;
    }
    if (family) {
        std::cout << "Beginning integration of family of solutions" << std::endl;
        double time1 = (double) clock() / CLOCKS_PER_SEC;
        TOV_Family sols(eos1, eos2, n_samples);
        double time2 = (double) clock() / CLOCKS_PER_SEC;
        std::cout << "Execution took " << time2-time1 << " seconds" << std::endl;
        std::cout << "DONE! Now printing a test lambda" << std::endl;
        std::cout << sols.lambda_from_mass(2.0) << std::endl;
    }
    else {
        TOV_result result = model.integrate_two_fluid_tov(e01, e02);
        model.print_result(result);
    }
    return 0;
}
