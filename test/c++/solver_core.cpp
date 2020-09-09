#include <triqs/utility/callbacks.hpp>
#include <iostream>

#include "moves/change_order.hpp"
#include "moves/change_order_alpha_sym.hpp"
#include "moves/change_order_cluster_groupings.hpp"
#include "moves/change_order_haar_scattering.hpp"
#include "moves/change_order_cluster_groupings_alpha_sym.hpp"
#include "moves/change_order_cluster_last_groupings_alpha_sym.hpp"
#include "moves/change_order_haar_scattering_alpha_sym.hpp"

#include "measures/measure_density_lo.hpp"
#include "measures/measure_density_mix.hpp"
#include "measures/measure_density_mix_groupings.hpp"
#include "measures/measure_histograms.hpp"
#include "measures/measure_correlation_lo.hpp"
#include "measures/measure_correlation_mix.hpp"
#include "measures/measure_times.hpp"

#include "solver_core.hpp"
#include "configuration.hpp"

using namespace triqs::gfs;
//namespace mpi = triqs::mpi;
using triqs::utility::mindex;

// ------------ The main class of the solver ------------------------

void solver_core::solve(solve_parameters_t const &params) {
 // Warning - Exception nb de points?

 // We choose the vectors to be always of size n + 1, even if the sampling can begin at min_perturbation_order
 d.resize(params.perturbation_order + 1);
 eta.resize(params.perturbation_order + 1);
 proba.resize(params.perturbation_order + 1);
 av_sign.resize(params.perturbation_order + 1);

 d() = 0;
 eta = d;
 proba = d;
 av_sign = d;

 // Needed containers
 auto g0_shift_R = params.g0_R;
 auto g0_shift_K = params.g0_K;
 g0_lesser_w = params.g0_R;
 g0_greater_w = g0_lesser_w;

 // Alpha shift + g>, g< definition
 for (auto w : g0_shift_R.mesh()) {
  g0_shift_R[w](0,0) = 1 / (1 / params.g0_R[w](0,0) - params.U * params.alpha);

  dcomplex fac1 = g0_shift_R[w](0,0) / params.g0_R[w](0,0);
  g0_shift_K[w](0,0) = fac1 * conj(fac1) * params.g0_K[w](0,0);

  dcomplex fac2 = 2_j * imag(g0_shift_R[w](0, 0));
  g0_lesser_w[w](0, 0) = (g0_shift_K[w](0,0) - fac2) / 2;
  g0_greater_w[w](0, 0) = (g0_shift_K[w](0,0) + fac2) / 2;
 }

 g0_lesser_t = make_gf_from_fourier(g0_lesser_w);
 g0_greater_t = make_gf_from_fourier(g0_greater_w);

 auto diag = make_gf_from_fourier(g0_lesser_w, true)(0)(0,0);

 // A useful bool
 bool is_alpha_sym = (real(g0_greater_t(0.5 * params.tmax)(0,0) - g0_lesser_t(0.5 * params.tmax)(0,0)) < 1.e-10);

 // Construct configuration
 configuration config(g0_lesser_t, g0_greater_t, params.alpha, params.tmax, diag, params.perturbation_order);

 // If order = 0, do note enter the MC loop
 if (params.perturbation_order == 0) {
  d(0) = imag(diag);
  eta(0) = abs(d(0));
  proba(0) = 1;
  av_sign(0) = 1;
  return;
 }

 if (world.rank() == 0) {
  std::cout << "Welcome to the LO CTDet solver" << std::endl;
  if (params.measure == "n_lo") std::cout << "You measure the density with the LO algorithm" <<std::endl;
  else if (params.measure == "n_lo_cluster_groupings") std::cout << "You measure the density with the LO algm, grouping according to cluster rules" << std::endl;
  else if (params.measure == "n_lo_cluster_last_groupings") std::cout << "You measure the density with the LO algm, grouping according to cluster rules and flipping last l" << std::endl;
  else if (params.measure == "n_lo_haar_scattering") std::cout << "You measure the density with the LO algm, grouping according to haar scattering perm vec" << std::endl;
  else if (params.measure == "n_mix") std::cout << "You measure the density with the mix algorithm" <<std::endl;
  else if (params.measure == "n_mix_cluster_groupings") std::cout << "You measure the density with the mix algm, grouping according to cluster rules" << std::endl;
  else if (params.measure == "n_mix_cluster_last_groupings") std::cout << "You measure the density with the mix algm, grouping according to cluster rules and flipping last l" << std::endl;
  else if (params.measure == "times") std::cout << "You measure the times visited by the MC" << std::endl;
  else TRIQS_RUNTIME_ERROR << "This is not a measurable quantity";
 }

 // Construct a MC - sign of first config ???
 triqs::mc_tools::mc_generic<std::complex<double>> MC(params.random_name, params.random_seed, params.verbosity);

 double sign_config = 1.0;

 if (params.min_perturbation_order > 1){
   auto v1 = get_random_lo_vertex(params.tmax, &(MC.get_rng()), params.cauchy_t, params.cauchy_a);
   auto v2 = get_random_lo_vertex(params.tmax, &(MC.get_rng()), params.cauchy_t, params.cauchy_a);
   
   // Prevent the sign to be zero
   v1.itau = 1;
   v2.itau = 1;

   if (v1.t < v2.t) {
     config.vertices = {v1,v2};
     config.times = {v1.t,v2.t};
   }
   else {
     config.vertices = {v2,v1};
     config.times = {v2.t, v1.t};
   }

   config.order = 2;
   
   std::complex<double> det = 1.0;
   for (int sigma = 0; sigma < 2; sigma++) det *= config.compute_det(sigma, config.vertices);
   config.det = det;

   config.sum_keldysh_indices = config.recompute_sum_keldysh_indices(config.times);
   
   std::complex<double> i{0,1};
   sign_config = (real(i * det) > 0) ? 1 : -1;
 }


 // Register moves - a bit complicated to read because the code is so modular
 if ((params.measure == "n_lo") || (params.measure == "n_mix") || (params.measure == "times")) {
   if (is_alpha_sym){
     MC.add_move(change_order_alpha_sym{&config, &(MC.get_rng()), params.perturbation_order, params.min_perturbation_order, params.tmax, params.lambda*params.U, params.cauchy_t,
                        params.cauchy_a, diag}, "Change order");
   }
   else {
     MC.add_move(change_order{&config, &(MC.get_rng()), params.perturbation_order, params.min_perturbation_order, params.tmax, params.lambda*params.U, params.cauchy_t,
                        params.cauchy_a, diag}, "Change order");
   }
 }

 else if ((params.measure == "n_lo_cluster_groupings") || (params.measure == "n_mix_cluster_groupings")) {
   if (is_alpha_sym) {
     MC.add_move(change_order_cluster_groupings_alpha_sym{&config, &(MC.get_rng()), params.perturbation_order, params.min_perturbation_order, params.tmax, params.lambda*params.U, params.cauchy_t,
                        params.cauchy_a, diag}, "Change order");
   }
   else {
     MC.add_move(change_order_cluster_groupings{&config, &(MC.get_rng()), params.perturbation_order, params.min_perturbation_order, params.tmax, params.lambda*params.U, params.cauchy_t,
                        params.cauchy_a, diag}, "Change order");
   }
 }

 else if ((params.measure == "n_lo_cluster_last_groupings") || (params.measure == "m_mix_cluster_last_groupings")) {
   if (is_alpha_sym) {
     MC.add_move(change_order_cluster_last_groupings_alpha_sym{&config, &(MC.get_rng()), params.perturbation_order, params.min_perturbation_order, params.tmax, params.lambda*params.U, params.cauchy_t,
                        params.cauchy_a, diag}, "Change order");
   }
   else {
     //MC.add_move(change_order_cluster_groupings{&config, &(MC.get_rng()), params.perturbation_order, params.min_perturbation_order, params.tmax, params.lambda*params.U, params.cauchy_t,
     //                   params.cauchy_a, diag}, "Change order");
   }
 }

 else if ((params.measure == "n_lo_haar_scattering")) {
   if (is_alpha_sym) {
     MC.add_move(change_order_haar_scattering_alpha_sym{&config, &(MC.get_rng()), params.perturbation_order, params.min_perturbation_order, params.tmax, params.lambda*params.U, params.cauchy_t,
                        params.cauchy_a, diag, params.perm_vecs}, "Change order");
   }
   else {
     MC.add_move(change_order_haar_scattering{&config, &(MC.get_rng()), params.perturbation_order, params.min_perturbation_order, params.tmax, params.lambda*params.U, params.cauchy_t,
                        params.cauchy_a, diag, params.perm_vecs}, "Change order");
   }
 }

 // Register measures
 // Note that we disable timer in all measures because system might spent its time calling the clock - updates being polynomial hence quick
 if ((params.measure == "n_lo") || (params.measure == "n_lo_cluster_groupings") || (params.measure == "n_lo_cluster_last_groupings") || (params.measure == "n_lo_haar_scattering")) {
   MC.add_measure(measure_density_lo{&config, &d, &eta, &proba, &av_sign, params.perturbation_order, params.lambda}, "Density measurement", false);
   if (params.measure_correlation) MC.add_measure(measure_correlation_lo{&config, correlation, params.perturbation_order,
                                                                    params.n_cycles, params.lambda}, "Opt - Correlation measurement", false);
 }
 else if (params.measure == "n_mix") {
   MC.add_measure(measure_density_mix{&config, &d, &eta, &proba, &av_sign, params.perturbation_order, params.tmax, params.lambda, is_alpha_sym}, "Density measurement", false);
   if (params.measure_correlation) MC.add_measure(measure_correlation_lo{&config, correlation, params.perturbation_order,
                                                                    params.n_cycles, params.lambda}, "Opt - Correlation measurement", false);
 }
 else if ((params.measure == "n_mix_cluster_groupings") || (params.measure == "n_mix_cluster_last_groupings")) {
   MC.add_measure(measure_density_mix_groupings{&config, &d, &eta, &proba, &av_sign, params.perturbation_order, params.tmax, params.lambda, is_alpha_sym}, "Density measurement", false);
   if (params.measure_correlation) MC.add_measure(measure_correlation_lo{&config, correlation, params.perturbation_order,
                                                                    params.n_cycles, params.lambda}, "Opt - Correlation measurement", false);
 }

 else if (params.measure == "times") MC.add_measure(measure_times{&config, times, params.n_times, params.perturbation_order}, "Times array measurement", false);

 if (params.measure_histos) MC.add_measure(measure_histograms{&config, times_histogram, weights_histogram, params.tmax,
                                                           params.perturbation_order, params.lambda}, "Opt - Histograms measurement", false); 
 
 // Run
 MC.warmup(params.n_warmup_cycles, params.length_cycle, triqs::utility::clock_callback(-1), sign_config);
 MC.accumulate(params.n_cycles, params.length_cycle, triqs::utility::clock_callback(params.max_MC_time));

 // Collect results
 MC.collect_results(world);
}
