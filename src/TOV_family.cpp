#include "../include/TOV_family.hpp"
#include <boost/math/interpolators/pchip.hpp>
#include <boost/math/tools/minima.hpp>
#include <boost/bind.hpp>
#include "../include/conversions.hpp"
#include "../include/spline.h"
#include <bits/stdc++.h>
#include <algorithm>


using boost::math::interpolators::pchip;
using boost::math::tools::brent_find_minima;
using boost::bind;

const int GRID_SIZE = 5;


template<typename T>
std::vector<double> linspace(T start_in, T end_in, int num_in)
{

  std::vector<double> linspaced;

  double start = static_cast<double>(start_in);
  double end = static_cast<double>(end_in);
  double num = static_cast<double>(num_in);

  if (num == 0) { return linspaced; }
  if (num == 1) 
    {
      linspaced.push_back(start);
      return linspaced;
    }

  double delta = (end - start) / (num - 1);

  for(int i=0; i < num-1; ++i)
    {
      linspaced.push_back(start + delta * i);
    }
  linspaced.push_back(end); // I want to ensure that start and end
                            // are exactly the same as the input
  return linspaced;
}

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

//inline void print_return()


TOV_Family::TOV_Family(EOS_Tabular& eos1, EOS_Poly& eos2, const int n_samples = 200) : eos1(eos1), eos2(eos2), n_samples(n_samples) {
    const double e1_max = eos1.e_max * CONVERSION::geom_to_dens_GeV_fm3;
    const double e1_min = eos1.e_min * CONVERSION::geom_to_dens_GeV_fm3;
    const double e2_max = e1_max*10;
    const double e2_min = e1_min*10;
    std::vector<double> e1s = linspace(e1_min, e1_max, n_samples);
    std::vector<double> e2s = linspace(e2_min, e2_max, n_samples);
    initialize_splines(e1s, e2s);
}

TOV_Family::TOV_Family(EOS_Tabular& eos1, EOS_Poly& eos2, std::string filename) : eos1(eos1), eos2(eos2), n_samples(n_samples) {
    // TODO: must use filesystem::path to deal with file IO
    TwoFluid_TOV model(eos1, eos2);
    const std::string table_path = "fx_tables/";
    const std::string full_filename = table_path+filename;
    std::ifstream f_eos;
    f_eos.open(full_filename);
    if (!f_eos.is_open()) {
        std::cout << "Cannot open file" << std::endl;
        exit(0);
    }
    f_eos.seekg(0);
    std::vector<double> e1s, e2s;
    while(true) {
        double e1, e2;
        f_eos >> e1;
        f_eos >> e2;
        if (f_eos.eof()) {break;}
        e1s.push_back(e1);
        e2s.push_back(e2);
    }
    initialize_splines(e1s, e2s);
}

TOV_Family::TOV_Family(EOS_Tabular& eos1, EOS_Poly& eos2, std::vector<double>& e1s, std::vector<double>& e2s) : eos1(eos1), eos2(eos2) {
    n_samples = e1s.size();
    initialize_splines(e1s, e2s);
}

TOV_Family::TOV_Family(EOS_Tabular& eos1, EOS_Poly& eos2) : eos1(eos1), eos2(eos2), n_samples(N_SAMPLE) {
    n_samples = N_SAMPLE;
    const double e1_min = 0.15;//eos1.e_min * CONVERSION::geom_to_dens_GeV_fm3;
    const double e1_max = eos1.e_max * CONVERSION::geom_to_dens_GeV_fm3;
    e1s = linspace(e1_min, e1_max, n_samples);
    
}

double TOV_Family::calc_F_chi(const double& K, const double& e1, const double& e2) {
    eos2.K = K;
    TwoFluid_TOV model(eos1, eos2);
    TOV_result res = model.integrate_two_fluid_tov(e1, e2);
    return res.F_chi;
}

double TOV_Family::calc_e2_from_F_chi(const double& K, const double& F_chi, const double& e1) {
    const int double_bits = std::numeric_limits<double>::digits;
    const double err_max = 1e-1;
    double lower_bound = 0.01;
    double upper_bound = 130;
    /*
    if (F_chi != 1) {
        double factor = pow(5.0, log10(K/1e7)) * (1e9/K) * (F_chi / (1.0-F_chi));
        if (e1 < 1) {
            upper_bound = factor * e1;
        }
        else {
            upper_bound = factor * pow(e1,2.0);
        }
    }
    else {
        upper_bound = (1e10/K) * pow(e1,2.0);
    }
    lower_bound = 1/3*upper_bound;
    */
    auto F_chi_minima = [&F_chi, &K, &e1, this] (double e2) {return 100*fabs(calc_F_chi(K, e1, e2)-F_chi);};
    std::pair<double, double> brent_root = brent_find_minima(F_chi_minima, lower_bound, upper_bound, double_bits);
    
    while (brent_root.second > err_max) {
        double try_e2 = brent_root.first;
        double err_e2 = 100*calc_F_chi(K, e1, try_e2)-F_chi;
        double center_point = try_e2;
        double dim_e2 = 0;
        double delta;
        if (try_e2 >= 1) {
            while (try_e2 > 10) {
                try_e2 /= 10;
                dim_e2 += 1;
            }
        }
        else {
            while(try_e2 < 1) {
                try_e2 *= 10;
                dim_e2 -= 1;
            }
        }
        delta = pow(10, dim_e2-1);
        if (err_e2 > 0) {
            //std::cout << "decreasing lower bound "; 
            lower_bound = center_point - 5*delta;
            upper_bound = center_point + delta;
        }
        else {
            //std::cout << "increasing upper bound ";
            lower_bound = center_point - delta;
            upper_bound = center_point + 5*delta;
        }
        /*
        std::cout << "Bounds: [" << lower_bound << ", " << upper_bound << "]";
        std::cout << " e1=" << e1;
        std::cout << " e2=" << center_point;
        std::cout << " F_chi=" << F_chi;
        std::cout << " delta=" << delta;
        std::cout << " dim_e2=" << dim_e2;
        std::cout << " err=" << err_e2 << "\r";
        std::cout.flush();
        */
        brent_root = brent_find_minima(F_chi_minima, lower_bound, upper_bound, double_bits);
    }
    return brent_root.first;
}

void TOV_Family::add_to_vector(double K, double F_chi, int index_start, int index_end) {
    for (int index=index_start; index < index_end; index++) {
        e2s[index] = calc_e2_from_F_chi(K, F_chi, e1s[index]);
    }
}

double TOV_Family::calc_lambda(const double K, const double F_chi, const double mass) {
    //std::vector<double> e1s = linspace(0.15, 3.0, N_SAMPLE);
    eos2.K = K;
    e2s = std::vector<double>(N_SAMPLE);
    for (int i = 0; i < N_SAMPLE; i++) {
        double e2 = calc_e2_from_F_chi(K, F_chi, e1s[i]);
        e2s[i] = e2;
    }
    initialize_splines(e1s, e2s);
    return lambda_from_mass(mass);
}

double TOV_Family::calc_lambda_parallel(const double K, const double F_chi, const double mass, const int n_threads) {
    //const int n_threads = 20;
    //std::cout << "Starting lambda calculation" << std::endl;
    if (N_SAMPLE % n_threads != 0) {
        std::cerr << "Number of threads must divide " << N_SAMPLE << std::endl;
        exit(1);
    }
    const int part = int(N_SAMPLE / n_threads);
    eos2.K = K;
    if (F_chi != 0) {
        e2s = std::vector<double>(N_SAMPLE);
        std::vector<std::thread> threads;
        int start = 0;
        for (int j = 0; j < n_threads; j++) {
            int index_start = start;
            int index_end = start+part;
            threads.push_back(std::thread(&TOV_Family::add_to_vector, this, K, F_chi, index_start, index_end));
            start = index_end;
            /*
            if (j >= 100) {
                threads[j-100].join();
            }
            */
        };
        for (auto& t : threads) {t.join();};
    }
    else {
        e2s = std::vector<double>(N_SAMPLE, 0);
    }
    //for (int i = 100; i < n_threads; i++) { threads[i].join();}
    //threads.clear();
    initialize_splines(e1s, e2s);
    double lambda = lambda_from_mass(mass);
    if (lambda > 5000) {
        //std::cout << "Wrong extrapolation for lambda: " << lambda << std::endl;
        //std::cout << "Setting lambda to max value (5000)" << std::endl;
        lambda = 5000;
    }
    else { 
        if (lambda < 0) {
            lambda = 0;
        }
        //std::cout << "Lambda: " << lambda << std::endl;
    }
    //std::cout << "Info: F_chi=" << F_chi << " K=" << K << " Mass=" << mass << std::endl;
    return lambda;
}

double TOV_Family::calc_lambda_normalized(double K_norm, double F_chi_norm, double mass_norm) {
    double K = pow(10, normalize_inverse(K_norm, log10(K_MIN), log10(K_MAX)));
    double F_chi = normalize_inverse(F_chi_norm, F_CHI_MIN, F_CHI_MAX);
    double mass = normalize_inverse(mass_norm, M_MIN, M_MAX);
    return calc_lambda(K, F_chi, mass);
}

void TOV_Family::generate_lambda_f_points() {
    std::vector<double> Ks = linspace(0.0, 1.0, GRID_SIZE);
    std::vector<double> F_chis = linspace(0.0, 1.0, GRID_SIZE);
    std::vector<double> Ms = linspace(0.0, 1.0, GRID_SIZE);
    Grid3D data_grid(GRID_SIZE);
    int index = 0;
    float progress = 0;
    std::cout << "Generating points" << std::endl;
    for (int i = 0; i < GRID_SIZE; i++) {
        for (int j = 0; j < GRID_SIZE; j++) {
            for (int k = 0; k < GRID_SIZE; k++) {
                double lambda_ijk = calc_lambda_normalized(Ks[k], F_chis[j], Ms[i]);
                data_grid.insert_point(index, {Ks[k], F_chis[j], Ms[i], lambda_ijk});
                reset_state();
                progress = index/pow(GRID_SIZE, 3);
                index++;
                std::cout << "Progress: ";
                std::cout << float(progress*100.0) << "%\r";
                std::cout.flush();
            }
        }
    };
    std::cout<<std::endl;
    data_grid.write_to_file("APR_POLY_N20.dat");

}

void TOV_Family::reset_state() {
    Rs.clear();
    Ms.clear();
    k2s.clear();
}

void TOV_Family::initialize_splines(const std::vector<double>& e1s, const std::vector<double>& e2s) {
    //std::cout << "Initializing splines" << std::endl;
    TwoFluid_TOV model(eos1, eos2);
    for (int i = 0; i < n_samples; i++) {
        TOV_result res = model.integrate_two_fluid_tov(e1s[i], e2s[i]);
        Rs.push_back(res.R);
        Ms.push_back(res.M);
        k2s.push_back(res.k2);
        model.reset_state();
    }
    sort_mass();
    remove_equal_entries();
    std::vector<double> Ms_copy = Ms;
    std::vector<double> Rs_copy = Rs;
    std::vector<double> k2s_copy = k2s;
    //tk::spline r_of_m(Ms, Rs);
    //tk::spline k2_of_m(Ms, k2s);
    r_m.setData(Ms, Rs);
    k2_m.setData(Ms, k2s);

    //r_m = pchip(std::move(Ms_copy), std::move(Rs_copy));
    Ms_copy = Ms;
    //k2_m = pchip(std::move(Ms_copy), std::move(k2s_copy));
    //r_m = r_of_m;
    //k2_m = k2_of_m;
    //std::cout << "Splines initialized" << std::endl;
}

double TOV_Family::radius_from_mass(double mass) {
    mass *= CONVERSION::Msun_to_mass_geom;
    return r_m(mass);
}

double TOV_Family::k2_from_mass(double mass) {
    mass *= CONVERSION::Msun_to_mass_geom;
    return k2_m(mass);
}

double TOV_Family::lambda_from_mass(double mass) {
    const double radius = radius_from_mass(mass);
    const double k2 = k2_from_mass(mass);
    mass *= CONVERSION::Msun_to_mass_geom;
    const double C = mass / radius;
    const double lambda = (2.0 / 3.0) * k2 / pow(C, 5);
    return lambda;
}

void TOV_Family::write_to_file() {
    if (Rs.size() > 0) {
        // TODO
        return;
    }
    else {
        std::cout << "There is no available data to write to file" << std::endl;
    }
}

void TOV_Family::sort_mass() {
    for (int j = 0; j < n_samples-1; j++) {
        for (int i = 0; i < n_samples-1-j; i++) {
            if (Ms[i+1] < Ms[i]) {
                std::swap(Ms[i], Ms[i+1]);
                std::swap(Rs[i], Rs[i+1]);
                std::swap(k2s[i], k2s[i+1]);
            }
        }
    }
}

void TOV_Family::remove_equal_entries() {
    //std::cout << "Removing" << std::endl;
    for (int i = n_samples-1; i >= 0; i--) {
        if(Ms[i] == Ms[i-1]) {
            //std::cout << Ms[i] << " " << Ms[i-1] << std::endl;

        }  
    }
    n_samples = Ms.size();
}