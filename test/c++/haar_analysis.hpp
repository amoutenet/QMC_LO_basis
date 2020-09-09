#pragma once
#include <triqs/gfs.hpp>
#include <triqs/arrays.hpp>
#include "bare_green_fn.hpp"
#include "lo_vertex.hpp"
#include "my_determinant.hpp"

using namespace triqs::gfs;
using namespace triqs::arrays;

std::pair<array<int, 1>, array<double, 1>>  haar_analysis(array<double, 1> times, double U, gf<retime> g0_lesser, gf<retime> g0_greater, dcomplex diag, dcomplex alpha, double t_max){
  
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
  array<double, 1> weights(n_weights / 2); // - (1 << order));
  array<int, 1> indices(weights.size());

  // Matrices for the dets
  matrix<std::complex<double>> A(order + 1, order + 1);

  int pos = 0;
  for (int it=0; it < n_weights; it++){
    // First case that we don't consider: the last itau is zeo
    if ((it >> 1) & 1 == 1) {

      // Alpha sym case
      /*bool b = true;
      int i = 0;
      while (b and (i < order - 1)) {
        b = (((it >> 2*i + 3) & 1) == ((order+1)%2));
        i += 1;
      }
      */

      //if (!b) {
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
  
        indices(pos) = it;
        weights(pos) = - std::pow(U, order) * real(coeffs[order %4] * _det) / (1 << (order + 1));
  
        pos += 1;
      //}
    }
  }

  return {indices, weights};
}
  
