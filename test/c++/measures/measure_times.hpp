#pragma once
#include <triqs/mc_tools.hpp>
#include <triqs/arrays.hpp>
#include <triqs/statistics/histograms.hpp>

using namespace triqs::gfs;
//namespace mpi = triqs::mpi;

struct measure_times {

 configuration* config; // Pointer to the MC qmc_data_t
 std::optional<triqs::arrays::matrix<double>>& times;
 int max_order;

 long it = 0;
 long n_times;

 measure_times(configuration* config_, std::optional<triqs::arrays::matrix<double>>& times_, long n_, int max_order_): 
   config(config_), max_order(max_order_), times(times_), n_times(n_) {
  std::cout << "Creation" << std::endl;
   times = triqs::arrays::matrix<double>(n_times, max_order);
   (*times)() = 0;
 }

 void accumulate(dcomplex sign) {

  if ((config->order == max_order) && (it < n_times)) {
    auto sli = (*times)(it, range());
    for (int i = 0; i < config->order; i++){
      sli(i) = config->vertices[i].t;
    }
   it++;
  }
 }

 void collect_results(mpi::communicator const &c) {
  *times = all_reduce(*times, c);
  if (it != n_times) TRIQS_RUNTIME_ERROR << "Problem in array size, it is " << it << " and n_times is " << n_times;
 }

};
