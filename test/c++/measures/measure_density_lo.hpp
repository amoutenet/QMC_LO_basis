#pragma once
#include <triqs/mc_tools.hpp>
#include <triqs/arrays.hpp>
#include <numeric>

using namespace triqs::gfs;

struct measure_density_lo {

 configuration* config; // Pointer to the MC qmc_data_t
 vector<double>& d;
 vector<double>& eta;
 vector<double>& proba;
 vector<double>& av_sign;

 int max_order;
 std::vector<double> lambda_coeffs;

 measure_density_lo(configuration* config_, vector<double>* d_, vector<double>* e_, vector<double>* p_, vector<double>* s_, int max_order_, double lambda_)
    : config(config_), d(*d_), eta(*e_), proba(*p_), av_sign(*s_), max_order(max_order_), lambda_coeffs(max_order_ + 1) {
  d() = 0;
  eta() = 0;
  proba() = 0;
  av_sign() = 0;

  lambda_coeffs[0]= 1;
  for (int i=1; i<lambda_coeffs.size(); i++) lambda_coeffs[i] = lambda_coeffs[i-1] / lambda_;
 }

 void accumulate(dcomplex sign) {
  d(config->order) += lambda_coeffs[config->order] * real(sign);
  eta(config->order) += lambda_coeffs[config->order];

  proba(config->order) += 1;
  av_sign(config->order) += real(sign);
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
  proba() /= n_measures;
 }
};
