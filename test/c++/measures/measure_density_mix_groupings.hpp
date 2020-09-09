#pragma once
#include <triqs/mc_tools.hpp>
#include <triqs/arrays.hpp>
#include <cmath>

using namespace triqs::gfs;
//namespace mpi = triqs::mpi;

struct measure_density_mix_groupings {

 configuration* config; // Pointer to the MC qmc_data_t
 vector<double>& d;
 vector<double>& eta;
 vector<double>& proba;
 vector<double>& av_sign;

 int max_order;
 double t_max;
 std::vector<double> lambda_coeffs;

 std::complex<double> i{0,1};
 std::vector<std::complex<double>> coeffs;

 std::vector<int> non_zero;

 measure_density_mix_groupings(configuration* config_, vector<double>* d_, vector<double>* e_, vector<double>* p_, vector<double>* s_, int max_order_, double tm_, double lambda_, bool is_alpha_sym)
    : config(config_), d(*d_), eta(*e_), proba(*p_), av_sign(*s_), max_order(max_order_), t_max(tm_), lambda_coeffs(max_order_ + 1), coeffs{i,{-1,0},-i,{1,0}}, non_zero(max_order_ + 1) {

  if (!is_alpha_sym) {
    non_zero[0] = 1;
       
    for (int i=1; i<non_zero.size(); i++) {
      int pow = 1 << i;
      non_zero[i] = pow * pow / 2;
    }
  } 
  else{
    non_zero[0] = 1;
    non_zero[1] = 2; // memory access ok - you should enter this fn only if config->order > 0
        
    for (int i=2; i<non_zero.size(); i++) {
      int pow = 1 << i;
      non_zero[i] = pow * pow / 2 - pow;
    }
  }

  d() = 0;
  eta() = 0;
  proba() = 0;
  av_sign() = 0;

  lambda_coeffs[0] = 1;
  for (int i=1; i<lambda_coeffs.size(); i++) lambda_coeffs[i] = lambda_coeffs[i-1] / lambda_;
 }

 void accumulate(dcomplex sign) {
  std::vector<double> V(config->vertices.size());
  for (int i = 0; i < V.size(); i++) V[i] = config->vertices[i].t;
  config->times = V;
  
  auto sum_dets = config->recompute_sum_keldysh_indices(config->times);
  config->sum_keldysh_indices = sum_dets;

  int k = config->order;
  double lo_weight = (-1) * real(coeffs[config->order %4] * config->det) / (1 << (k+1));

  auto x = std::isfinite(real(sum_dets)) ? real(sum_dets) : 0;

  d(k) += lambda_coeffs[k] * x / (std::abs(lo_weight) * non_zero[k]);
  eta(k) += lambda_coeffs[k];
  proba(k) += 1;
  av_sign(k) += real(sign);
  }

 void collect_results(mpi::communicator c) {
  d = all_reduce(d, c);
  eta = all_reduce(eta, c);
  proba = all_reduce(proba, c);
  av_sign = all_reduce(av_sign, c);

  double n_measures = std::accumulate(proba.begin(), proba.end(), double{0});
  d() /= n_measures;
  eta() /= n_measures;

  for (int i=0; i<av_sign.size(); i++) av_sign(i) /= proba(i);
  proba /= n_measures;
 }
};
