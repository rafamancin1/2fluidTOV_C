#include "../include/TOV_family.hpp"
//#include <boost/math/interpolators/pchip.hpp>
#include <boost/math/tools/minima.hpp>
#include <boost/bind.hpp>
#include "../include/conversions.hpp"
#include <bits/stdc++.h>
#include <algorithm>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_min.h>
#include <gsl/gsl_roots.h>
#include <boost/math/tools/roots.hpp> 

//using boost::math::interpolators::pchip;
using boost::math::tools::brent_find_minima;
using boost::bind;
EOS_Poly null_eos("Polytropic", 0.0, 0.0);

const int GRID_SIZE = 30;

template< typename F >  class gsl_function_pp : public gsl_function {
 public:
 gsl_function_pp(const F& func) : _func(func) {
   function = &gsl_function_pp::invoke;
   params=this;
 }
 private:
 const F& _func;
 static double invoke(double x, void *params) {
 return static_cast<gsl_function_pp*>(params)->_func(x);
 }
 };


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


TOV_Family::TOV_Family(EOS_Tabular& eos1, EOS_Analytic& eos2, int n_samples)
    : eos1(eos1), eos2(eos2), n_samples(n_samples) {
    const double e1_max = eos1.e_max * CONVERSION::geom_to_dens_GeV_fm3;
    const double e1_min = eos1.e_min * CONVERSION::geom_to_dens_GeV_fm3;
    const double e2_max = e1_max*10;
    const double e2_min = e1_min*10;
    std::vector<double> e1s = linspace(e1_min, e1_max, n_samples);
    std::vector<double> e2s = linspace(e2_min, e2_max, n_samples);
    initialize_splines(e1s, e2s);
}

TOV_Family::TOV_Family(EOS_Tabular& eos1, EOS_Analytic& eos2, std::string filename)
    : eos1(eos1), eos2(eos2), n_samples(n_samples) {
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

TOV_Family::TOV_Family(EOS_Tabular& eos1, EOS_Analytic& eos2, std::vector<double>& e1s, std::vector<double>& e2s)
    : eos1(eos1), eos2(eos2) {
    n_samples = e1s.size();
    initialize_splines(e1s, e2s);
}

TOV_Family::TOV_Family(EOS_Tabular& eos1, EOS_Analytic& eos2)
    : eos1(eos1), eos2(eos2), n_samples(N_SAMPLE) {
    n_samples = N_SAMPLE;
    const double e1_min = 0.55;//eos1.e_min * CONVERSION::geom_to_dens_GeV_fm3;
    const double e1_max = 2.4;
    //std::cout << "Initializing TOV Family with e_min=" << e1_min << " and e_max=" << e1_max << std::endl;
    assert(e1_max > 0 && e1_max > e1_min);
    //e1s = linspace(e1_min, e1_max, n_samples);
}

TOV_Family::TOV_Family(EOS_Tabular& eos1, EOS_Analytic& eos2, double F_chi)
    : eos1(eos1), eos2(eos2), n_samples(N_SAMPLE) {
    n_samples = N_SAMPLE;
    const double e1_min = 0.15;
    const double e1_max = 2.4; //eos1.e_max * CONVERSION::geom_to_fm3 - 0.2;
    //std::cout << e1_max << std::endl;
    assert(e1_max > 0);
    e1s = linspace(e1_min, e1_max, n_samples);
    generate_e2s_from_F_chi(F_chi);
    initialize_splines(e1s, e2s);
}

TOV_Family::TOV_Family(EOS_Tabular& eos1) : eos1(eos1), eos2(null_eos), n_samples(N_SAMPLE) {
    const double e1_min = 0.3;
    const double e1_max = 1.44;//eos1.e_max * CONVERSION::geom_to_fm3 - 0.1;
    assert(e1_max > 0);
    e1s = linspace(e1_min, e1_max, n_samples);
    e2s = std::vector<double>(n_samples, 0.0);
    initialize_splines(e1s, e2s);
}

double TOV_Family::calc_F_chi(const double& e1, const double& e2) {
    //eos2.K = K;
    TwoFluid_TOV model(eos1, eos2);
    TOV_result res = model.integrate_two_fluid_tov(e1, e2);
    return res.F_chi;
}

std::pair<double, double> TOV_Family::calc_e2_from_F_chi(const double& F_chi, const double& e1) {
    const int double_bits = std::numeric_limits<double>::digits;
    const double err_max = 1;
    double lower_bound = 0.0;
    double upper_bound = 300.0;
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
    auto F_chi_minima = [&F_chi, &e1, this] (double e2) {return 100*fabs(calc_F_chi(e1, e2)-F_chi);};
    std::pair<double, double> brent_root = brent_find_minima(F_chi_minima, lower_bound, upper_bound, double_bits);
    
    //if (brent_root.second > err_max) {std::cout << brent_root.second << std::endl;}; 
    return brent_root;
}

double TOV_Family::calc_e2_from_F_chi_v2(const double& F_chi, const double& e1) {
    if (F_chi == 0.0) {
        return 0.0;
    }
    const double err_max = 0.0;
    const double err_abs = 1e-6;
    double lower_bound = 0.0;
    double upper_bound = 10.0;
    int status;
    //double guess = lower_bound + (3 - sqrt(5)/2)*(upper_bound-lower_bound);
    double root;
    auto F_chi_minima = [&F_chi, &e1, this] (double e2) {return 100*(calc_F_chi(e1, e2)-F_chi);};
    gsl_function_pp<decltype(F_chi_minima)> Fp(F_chi_minima);
    const gsl_root_fsolver_type *T;

    gsl_root_fsolver *s;
    gsl_function *F = static_cast<gsl_function*>(&Fp);
    T = gsl_root_fsolver_brent;
    s = gsl_root_fsolver_alloc(T);
    gsl_error_handler_t *default_handler = gsl_set_error_handler_off();
    int max_tries = 100;
    int n_tries = 0;
    do { 
        status = gsl_root_fsolver_set(s, F, lower_bound, upper_bound);
        if (status == GSL_EINVAL) { // endpoints do not straddle means supplied bounds are invalid
            upper_bound *= 2.0;
        }
        n_tries++;
    } while (status == GSL_EINVAL && n_tries <= max_tries);
    if (status == GSL_EINVAL) {
        std::cout << "Could not find working upper bound."; 
    }
    do {
        status = gsl_root_fsolver_iterate(s);
        root = gsl_root_fsolver_root(s);
        upper_bound = gsl_root_fsolver_x_upper(s);
        lower_bound = gsl_root_fsolver_x_lower(s);
        
        status = gsl_min_test_interval(lower_bound, upper_bound, err_abs, err_max);
        
    } while (status == GSL_CONTINUE);
    if (status) {
        std::cout << "Status: " << status << std::endl;
        std::cout << "Root finding did not converge" << std::endl;
    }
    //std::cout << "err " << F_chi_minima(root) << std::endl;
    //std::cout << "err(est)" << upper_bound - lower_bound << std::endl;
    //std::cout << upper_bound << std::endl;
    //std::cout << root << std::endl;
    //std::cout << "err " << F_chi_minima(root) << std::endl;
    gsl_set_error_handler(default_handler); // to reset the error handler
    gsl_root_fsolver_free(s);
    return root;
}

void TOV_Family::add_to_vector(double F_chi, int index_start, int index_end) {
    for (int index=index_start; index < index_end; index++) {
        double root = calc_e2_from_F_chi_v2(F_chi, e1s[index]);
        e2s[index] = root;
        /*
        if (root.second > 1) {
            bad_indexes.push_back(index);
        }
        */
    }
}

void TOV_Family::generate_e2s_from_F_chi(double F_chi) {
    if (F_chi == 0.0) {
        e2s = std::vector<double>(n_samples, 0.0);
    }
    else {
        e2s = std::vector<double>(N_SAMPLE);
        for (int i = 0; i < N_SAMPLE; i++) {
            double e2 = calc_e2_from_F_chi_v2(F_chi, e1s[i]);
            e2s[i] = e2;
        }
    }
}

void TOV_Family::generate_e2s_from_F_chi_v2(double K, double F_chi) {
    e2s = std::vector<double>(N_SAMPLE);
    const double C = 11;
    for (int i = 0; i < N_SAMPLE; i++) {
        double e2 = C * F_chi * e1s[i]; //calc_e2_from_F_chi_v2(K, F_chi, e1s[i]);
        e2s[i] = e2;
    }
}

TOV_result TOV_Family::calc_lambda_and_mass_directly(double e1, double F_chi) {
    TwoFluid_TOV model(eos1, eos2);
    double e2 = calc_e2_from_F_chi_v2(F_chi, e1);
    TOV_result result = model.integrate_two_fluid_tov(e1, e2);
    return result;
}

TOV_result TOV_Family::calc_lambda_and_mass_directly_v2(double e1, double F_chi) {
    //eos2.K = K;
    TwoFluid_TOV model(eos1, eos2);
    double e2 = F_chi * e1; //calc_e2_from_F_chi_v2(eos2.K, F_chi, e1);
    TOV_result result = model.integrate_two_fluid_tov(e1, e2);
    return result;
}

double TOV_Family::calc_lambda(const double F_chi, const double mass) {
    //std::vector<double> e1s = linspace(0.15, 3.0, N_SAMPLE);
    //eos2.K = K;
    generate_e2s_from_F_chi(F_chi);
    initialize_splines(e1s, e2s);
    return lambda_from_mass(mass);
}

double TOV_Family::calc_lambda_parallel(const double F_chi, const double mass, const int n_threads) {
    //const int n_threads = 20;
    //std::cout << "Starting lambda calculation" << std::endl;
    if (N_SAMPLE % n_threads != 0) {
        std::cerr << "Number of threads must divide " << N_SAMPLE << std::endl;
        exit(1);
    }
    const int part = int(N_SAMPLE / n_threads);
    if (F_chi != 0) {
        e2s = std::vector<double>(N_SAMPLE);
        std::vector<std::thread> threads;
        int start = 0;
        for (int j = 0; j < n_threads; j++) {
            int index_start = start;
            int index_end = start+part;
            threads.push_back(std::thread(&TOV_Family::add_to_vector, this, F_chi, index_start, index_end));
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
    //for (int x : bad_indexes) {e1s.erase(e1s.begin()+x); e2s.erase(e2s.begin()+x); n_samples -= 1;}
    initialize_splines(e1s, e2s);
    double lambda = lambda_from_mass(mass);
    if (lambda > 5000) {
        //std::cout << "Wrong extrapolation for lambda: " << lambda << std::endl;
        //std::cout << "Setting lambda to max value (5000)" << std::endl;
        lambda = 5000;
    }
    else { 
        if (lambda < 0 || lambda != lambda) {
            lambda = 0;
        }
        //std::cout << "Lambda: " << lambda << std::endl;
    }
    //std::cout << "Info: F_chi=" << F_chi << " K=" << K << " Mass=" << mass << std::endl;
    return lambda;
}

double TOV_Family::calc_lambda_normalized(double K_norm, double F_chi_norm, double mass_norm) {
    const int N_THREADS = 10;
    double K = pow(10, normalize_inverse(K_norm, log10(K_MIN), log10(K_MAX)));
    double F_chi = normalize_inverse(F_chi_norm, F_CHI_MIN, F_CHI_MAX);
    double mass = normalize_inverse(mass_norm, M_MIN, M_MAX);
    double lambda = calc_lambda_parallel(F_chi, mass, N_THREADS);
    double lambda_norm = normalize(lambda, LAMBDA_MIN, LAMBDA_MAX); 

    return lambda_norm;
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
                std::cout << float(progress*100.0) << "% \r";
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
        RBs.push_back(res.R_B);
        RDs.push_back(res.R_D);
        lambdas.push_back(res.lambda);
        //std::cout << "E ";
        //std::cout << e1s[i] << std::endl;
        //std::cout << "Mass ";
        //std::cout << res.M * CONVERSION::mass_geom_to_Msun << std::endl;
        model.reset_state();
    }
    //sort_mass();
    //remove_equal_entries();
    r_m.setData(Ms, Rs);
    k2_m.setData(Ms, k2s);
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
    // TODO: try interpolating lambda directly
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

void TOV_Family::set_analytic_eos(EOS_Analytic& eos2) {
    this->eos2 = eos2;
}