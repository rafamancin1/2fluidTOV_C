#ifndef UTILS_HPP
#define UTILS_HPP
#include "TOV_family.hpp"
#include "../boost/math/tools/minima.hpp"
#include "../boost/bind.hpp"

using boost::math::tools::brent_find_minima;
using boost::bind;

const int N_SAMPLE = 200;
const int GRID_SIZE = 3;
const double M_MIN = 0.2;
const double M_MAX = 2.1;
const double K_MIN = 1e6;
const double K_MAX = 1e8;
const double F_chi_MIN = 0.01;
const double F_chi_MAX = 0.5;
const auto normalize = [](const double X, const double X_min, const double X_max) {return (X-X_min)/(X_max-X_min);};
const auto normalize_inverse = [](const double X_norm, const double X_min, const double X_max) {return (X_max-X_min)*X_norm + X_min;};

EOS_Tabular eos1("LALSimNeutronStarEOS_APR");

class Grid3D {

    std::vector<double> X;
    std::vector<double> Y;
    std::vector<double> Z;
    std::vector<double> W;
    int size;
    int total_size;

    public:
        Grid3D(int size) : size(size) {
            total_size = pow(size, 3);
            X = std::vector<double>(total_size);
            Y = std::vector<double>(total_size);
            Z = std::vector<double>(total_size);
            W = std::vector<double>(total_size);
        }
    
        void insert_point(const int index, const std::vector<double> point) {
            if (point.size() != 4) {
                std::cerr << "Point must be four-dimensional" << std::endl;
                exit(1);
            }
            else {
                X[index] = point[0];
                Y[index] = point[1];
                Z[index] = point[2];
                W[index] = point[3];
            }
        }
        std::vector<double> get_point(int index) {
            return {X[index], Y[index], Z[index], W[index]};
        }

        std::vector<double> get_point(int i, int j, int k) {
            int index = i*pow(size, 2) + j*size + k;
            return get_point(index);
        }

        void write_to_file(std::string filename) {
            std::ofstream file(filename);
            if (!file.is_open()) {
                std::cerr << "Could not open file" << std::endl;
                exit(1);
            }
            for (int i = 0; i < total_size; i++) {
                std::cout << "writing point (" << X[i] << " " << Y[i] << " " << Z[i] << " " << W[i] << ")" << std::endl;
                file << X[i] << " ";
                file << Y[i] << " ";
                file << Z[i] << " ";
                file << W[i] << " ";
                file << "\n";
            }
            file.close();
        }

};
/*
typedef struct Grid4D {
    std::vector<double> X;
    std::vector<double> Y;
    std::vector<double> Z;
    std::vector<double> W;
} Grid4D;
*/
double calc_F_chi(double K, double e1, double e2) {

    EOS_Poly eos2("Polytropic", K, 2);
    TwoFluid_TOV model(eos1, eos2);
    TOV_result res = model.integrate_two_fluid_tov(e1, e2);
    return res.F_chi;
}

double calc_e2_from_F_chi(double K, double F_chi, double e1) {
    const int double_bits = std::numeric_limits<double>::digits;
    double lower_bound;
    double upper_bound;
    // Setting bounds
    if (F_chi != 1) {
        double factor = pow(5,log10(K / 1e7)) * (1e9 / K) * (F_chi / (1-F_chi));
        if (e1 < 1) {
            upper_bound = factor * e1;
        }
        else {
            upper_bound = factor * pow(e1,2);
        }
    }
    else {
        upper_bound = (1e10 / K) * pow(e1,2);
    }
    lower_bound = 1/3 * upper_bound;
    auto F_chi_minima = [&F_chi, &K, &e1](double e2) {return fabs(calc_F_chi(K, e1, e2) - F_chi);};
    std::pair<double, double> brent_root = brent_find_minima(F_chi_minima, lower_bound, upper_bound, double_bits);
    return brent_root.first;
}

void test_calc_e2(double K, double F_chi, double e1) {
    double e2 = calc_e2_from_F_chi(K, F_chi, e1);
    double F_chi_right = calc_F_chi(K, e1, e2);
    std::cout << "F_chi target " << F_chi << std::endl;
    std::cout << "F_chi calculated from e2 found: " << F_chi_right << std::endl; 
    std::cout << fabs(F_chi-F_chi_right) << std::endl; 
}

double calc_lambda(double K, double F_chi, double mass) {
    std::vector<double> e1s = linspace(0.15, 3.0, N_SAMPLE);
    EOS_Tabular eos1("LALSimNeutronStarEOS_APR");
    EOS_Poly eos2("DM_Polytropic", K, 2);
    std::vector<double> e2s(200);
    for (int i = 0; i < N_SAMPLE; i++) {
        double e2 = calc_e2_from_F_chi(K, F_chi, e1s[i]);
        e2s[i] = e2;
    }
    TOV_Family sols(eos1, eos2, e1s, e2s);
    return sols.lambda_from_mass(mass);
}

double calc_lambda_normalized(double K_norm, double F_chi_norm, double mass_norm) {
    double K = pow(10, normalize_inverse(K_norm, log10(K_MIN), log10(K_MAX)));
    double F_chi = normalize_inverse(F_chi_norm, F_chi_MIN, F_chi_MAX);
    double mass = normalize_inverse(mass_norm, M_MIN, M_MAX);
    return calc_lambda(K, F_chi, mass);

}

void generate_lambda_f_points () {
    std::vector<double> Ks = linspace(0.0, 1.0, GRID_SIZE);
    std::vector<double> F_chis = linspace(0.0, 1.0, GRID_SIZE);
    std::vector<double> Ms = linspace(0.0, 1.0, GRID_SIZE);
    Grid3D data_grid(GRID_SIZE);
    int index = 0;
    std::cout << "Generating of points" << std::endl;
    for (int i = 0; i < GRID_SIZE; i++) {
        for (int j = 0; j < GRID_SIZE; j++) {
            for (int k = 0; k < GRID_SIZE; k++) {
                double lambda_ijk = calc_lambda_normalized(Ks[k], F_chis[j], Ms[i]);
                data_grid.insert_point(index, {Ks[k], F_chis[j], Ms[i], lambda_ijk});
                index++;
            }
        }
    };
    data_grid.write_to_file("APR_Poly_N20.dat");
}

#endif