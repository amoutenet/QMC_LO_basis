#pragma once
#include <triqs/mc_tools.hpp>
#include <triqs/arrays.hpp>
#include <triqs/statistics/histograms.hpp>

using namespace triqs::gfs;
//namespace mpi = triqs::mpi;

struct measure_histograms {

 configuration* config; // Pointer to the MC qmc_data_t
 //triqs::statistics::histogram* histo_perturbation_order;
 std::optional<triqs::statistics::histogram>& histo_times;
 std::optional<triqs::statistics::histogram>& histo_weights;

 int max_order;
 int power_two;
 double power_lambda;

 std::complex<double> i{0,1};
 std::vector<std::complex<double>> coeffs;

  measure_histograms(configuration* config_, std::optional<triqs::statistics::histogram>& hist_t, std::optional<triqs::statistics::histogram>& hist_w,
                     double tmax, int max_order_, double lambda_):
    config(config_), max_order(max_order_), power_two(1<<(max_order_+1)), power_lambda(std::pow(lambda_, max_order_)),
    histo_times(hist_t), histo_weights(hist_w), coeffs{i,{-1,0},-i,{1,0}} {
    //*histo_perturbation_order = {0, 1000};
   histo_times = {0, tmax, 1000};
   histo_weights = {-3, 3, 1000};
  }

  void accumulate(dcomplex sign) {
   if (config->order == max_order) {
    *histo_weights << (-1) * real(coeffs[config->order %4] * config->det) * power_lambda / power_two;

    for (int i = 0; i < config->order; ++i) {
     *histo_times << config->vertices[i].t;
    }
   }
  }

  void collect_results(mpi::communicator const &c) {

    *histo_times = all_reduce(*histo_times, c);
   *histo_weights = all_reduce(*histo_weights, c);

  }

};
