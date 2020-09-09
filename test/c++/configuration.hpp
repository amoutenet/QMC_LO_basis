#pragma once
#include <triqs/gfs.hpp>
#include "./bare_green_fn.hpp"
#include "./michel_determinant.hpp"
#include "./det_manip_corentin.hpp"
#include "./det_manip_triqs_v2.hpp"

/// Definition of a configuration

using namespace triqs::gfs;

struct configuration{
 int order =0; //current perturbation order

 dcomplex det;
 std::vector<lo_vertex> vertices;
 std::vector<double> times;

 dcomplex sum_keldysh_indices;
 dcomplex diag;

 michel_determinant detclass;
 lo_vertex vmax;
 double tmax;
 
 std::vector<bare_green_fn_lo> g0_lo;
 std::vector<bare_green_fn_pm> g0_pm;

 // We allocate ONCE the matrix A in memory and we'll consider views of it
 triqs::arrays::array<dcomplex, 2> A;
 std::vector<int> ri;

 configuration(gf<retime> const& g0_lesser, gf<retime> const& g0_greater, dcomplex alpha, double t_, dcomplex d_, int order): det(2*d_), A(order + 1, order + 1), ri(order + 1), diag(d_), tmax(t_) {

   // Define proxy vector
   for (int s = 0; s < 2; ++s) g0_lo.emplace_back(s, g0_lesser, g0_greater, alpha, t_, d_);
   for (int s = 0; s < 2; ++s) g0_pm.emplace_back(s, g0_lesser, g0_greater, alpha, t_, d_);
   sum_keldysh_indices = imag(d_);
   vmax = lo_vertex{t_,0,0};
 }

 dcomplex compute_det(int sigma, std::vector<lo_vertex> V);
 dcomplex compute_det(int sigma, std::vector<keldysh_pt> V);

 dcomplex recompute_sum_keldysh_indices(std::vector<double> V);

 std::pair<int, bool> get_flip_idx(std::vector<lo_vertex> V);
};
