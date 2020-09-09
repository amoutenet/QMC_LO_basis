#pragma once
#include <triqs/gfs.hpp>
#include <triqs/arrays.hpp>
#include "bare_green_fn.hpp"
#include "lo_vertex.hpp"
#include "my_determinant.hpp"

using namespace triqs::gfs;
using namespace triqs::arrays;

std::pair<array<int, 1>, array<double, 1>>  weight_analysis(array<double, 1> times, double U, gf<retime> g0_lesser, gf<retime> g0_greater, dcomplex diag, dcomplex alpha, double t_max){
  
  // Bare green's function proxy, a vector of size 2
  std::vector<bare_green_fn_lo> g0 {bare_green_fn_lo{0, g0_lesser, g0_greater, alpha, t_max, diag}, bare_green_fn_lo{1, g0_lesser, g0_greater, alpha, t_max, diag}};

  // Useful coefficients vector
  std::complex<double> i{0,1};
  std::vector<std::complex<double>> coeffs{i,{-1,0},-i,{1,0}};

  int order = times.size();
  array<lo_vertex, 1> vertices(order);
  my_determinant det(order);

  int n_weights = (1 << order) * (1 << order);

  // Construct the vector of (indices, weights) we're going to sort
  std::vector<std::pair<int, double>> weights(n_weights); 
 
  // Matrices for the dets
  matrix<std::complex<double>> A(order + 1, order + 1);

  for (int it=0; it < n_weights; it++){

    for (int n=0; n < order; n++) {
      vertices(n) = lo_vertex{times(order-1 - n), (it >> (2*n+1))&1, (it >> 2*n)&1};
    }

    auto v_max = lo_vertex{t_max, 0, 0};
    dcomplex _det = 1;

    // Compute the det
    for (int s = 0; s < 2; s++){

      for (int j = 0; j < order; j++){
	for (int k = 0; k < order; k++) A(k,j) = g0[s](vertices(k), vertices(j));
        if (s == 0){
	  A(order,j) = g0[s](v_max,vertices(j));
	  A(j,order) = g0[s](vertices(j),v_max);
	}
      }

      if (s == 0) A(order, order) = g0[s](v_max, v_max);

      _det *= (s == 0) ? det.det(&A, order+1) : det.det(&A, order);
    }

    weights[it] = {it, - std::pow(U, order) * real(coeffs[order %4] * _det) / (1 << (order + 1))};
  }

  std::sort(weights.begin(), weights.end(), [](auto const& x, auto const& y) {return std::abs(x.second) < std::abs(y.second);});

  array<int, 1> sorted_indices(n_weights);
  array<double, 1> sorted_weights(n_weights);

  for (int it=0; it < n_weights; it++) std::tie(sorted_indices(it), sorted_weights(it)) = weights[it];

  return {sorted_indices, sorted_weights};
}
  
