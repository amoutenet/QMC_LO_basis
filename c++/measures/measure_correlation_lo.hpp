#pragma once
#include <triqs/mc_tools.hpp>
#include <triqs/arrays.hpp>

using namespace triqs::gfs;
//namespace mpi = triqs::mpi;

struct measure_correlation_lo {

 configuration* config; // Pointer to the MC qmc_data_t
 std::optional<triqs::arrays::vector<double>>& correlation;
 int max_order;

 int power_two;
 double power_lambda;

 std::complex<double> i{0,1};
 std::vector<std::complex<double>> coeffs;

 long it = 0;

 measure_correlation_lo(configuration* config_, std::optional<triqs::arrays::vector<double>>& corr, int max_order_, long n_cycles, double lambda_)
   : config(config_), correlation(corr), max_order(max_order_), power_two(1<<(max_order_+1)),
     power_lambda(std::pow(lambda_, max_order_)), coeffs{i,{-1,0},-i,{1,0}} {
  correlation = triqs::arrays::vector<double>(n_cycles);
  
  (*correlation)() = 0;
 }

 void accumulate(dcomplex sign) {
  int k = config->order;

  if (k == max_order) {
   (*correlation)(it) = (-1) * real(coeffs[config->order %4] * config->det) * power_lambda / power_two;
   it++;
  }
 }

 void collect_results(mpi::communicator c) {
  std::cout << "Length ot the array " << it << std::endl;
  *correlation = all_reduce(*correlation, c);
 }
};

