#pragma once

#include <triqs/mc_tools.hpp>
#include <triqs/arrays.hpp>

using namespace triqs::gfs;
//namespace mpi = triqs::mpi;

struct measure_correlation_mix {

 configuration* config; // Pointer to the MC qmc_data_t
 std::optional<triqs::arrays::vector<double>>& correlation;
 int max_order;
 
 double power_lambda;
 long it = 0;

 measure_correlation_mix(configuration* config_, std::optional<triqs::arrays::vector<double>>& corr, int max_order_, long n_cycles, double lambda_)
   : config(config_), correlation(corr), max_order(max_order_), power_lambda(std::pow(lambda_, max_order_)) {
  
  correlation = triqs::arrays::vector<double>(n_cycles);
  (*correlation)() = 0;
 }

 void accumulate(dcomplex sign) {
  if (config->order == max_order) {
   (*correlation)(it) = real(config->sum_keldysh_indices) * power_lambda;
   it++;
  }
 }

 void collect_results(mpi::communicator c) {
  std::cout << "Length ot the array " << it << std::endl;
  *correlation = all_reduce(*correlation, c);
 }
};
