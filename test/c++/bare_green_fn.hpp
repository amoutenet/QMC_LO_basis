#pragma once
#include "lo_vertex.hpp"
#include "keldysh_pt.hpp"
#include <triqs/gfs.hpp>

using namespace triqs::gfs;

// --------------   G0 adaptor   --------------------------------------

/*
 * It is the function that appears in the calculation of the determinant (det_manip, cf below).
 */

struct bare_green_fn_lo {
 gf<retime> g0_lesser;
 gf<retime> g0_greater;

 int sigma; //spin
 gf<retime, matrix_valued> G_matrix;

 dcomplex alpha;
 double t_max;
 dcomplex diag;

 bare_green_fn_lo(int sigma_, gf<retime> const &g0_lesser_, gf<retime> const &g0_greater_, dcomplex alpha_, double t_, dcomplex d_) :
   g0_lesser(g0_lesser_), g0_greater(g0_greater_), sigma(sigma_), G_matrix{g0_lesser_.mesh(), {3, 3}}, alpha(alpha_), t_max(t_), diag(d_) {
  auto K = g0_lesser + g0_greater;

  auto pp = g0_lesser;
  for (auto const &t: pp.mesh()) {
   if (double(t) > 0) pp[t] = g0_greater[t];
  }

  auto R = pp - g0_lesser;
  auto A = pp - g0_greater;

  int eps = 1 - 2*sigma;

  //FIXME Olivier
  slice_target(G_matrix, range(0,1), range(0,1)) = R;
  slice_target(G_matrix, range(0,1), range(1,2)) = K;
  slice_target(G_matrix, range(0,1), range(2,3)) = R + eps * K;
  slice_target(G_matrix, range(1,2), range(0,1))() = 0;
  slice_target(G_matrix, range(1,2), range(1,2)) = A;
  slice_target(G_matrix, range(1,2), range(2,3)) = eps * A;
  slice_target(G_matrix, range(2,3), range(0,1)) = eps * R;
  slice_target(G_matrix, range(2,3), range(1,2)) = eps * K + A;
  slice_target(G_matrix, range(2,3), range(2,3)) = eps * (R + A) + K;
 }

 std::complex<double> operator()(lo_vertex v1, lo_vertex v2) const {
  if (v1.t == v2.t) {
   if ((v1.itau == sigma) and (v2.itau == sigma)) return 2*diag- ((v1.t == t_max) ? 0_j : 2_j * alpha);
   return 0;
  }

  //return G_matrix(v1.t - v2.t)(v1.itau == sigma ? 2 : v1.l, v2.itau == sigma ? 2 : v2.l);
  return G_matrix[triqs::gfs::closest_mesh_pt(v1.t - v2.t)](v1.itau == sigma ? 2 : v1.l, v2.itau == sigma ? 2 : v2.l);
 }

};

struct bare_green_fn_pm {

 gf<retime> g0_lesser;
 gf<retime> g0_greater;
 int sigma;

 dcomplex alpha;
 double t_max;
 dcomplex diag;

 bare_green_fn_pm(int sigma_, gf<retime> const &g0_lesser_, gf<retime> const &g0_greater_, dcomplex alpha_, double t_, dcomplex d_) :
   g0_lesser(g0_lesser_), g0_greater(g0_greater_), sigma(sigma_), alpha(alpha_), t_max(t_), diag(d_) {}

 dcomplex operator()(keldysh_pt const &a, keldysh_pt const &b) const {
  // at equal time, discard Keldysh index and use g_lesser
  // do not put alpha for the time_max even at equal time
  if (a.t == b.t) return diag - ((b.t == t_max) ? 0_j : 1_j * alpha);

  //  // mapping: is it lesser or greater?
  //  //  a    b    (a.time > b.time)   L/G ?
  //  //  0    0           1             G
  //  //  0    0           0             L
  //  //  1    1           1             L
  //  //  1    1           0             G
  //  //
  //  //  0    1           *             L
  //  //  1    0           *             G
  bool is_greater = (a.k == b.k ? (a.k xor (a.t > b.t)) : a.k);
  //return (is_greater ? g0_greater[triqs::gfs::closest_mesh_pt(a.t - b.t)](0,0) : g0_lesser[triqs::gfs::closest_mesh_pt(a.t - b.t)](0,0));
  return (is_greater ? g0_greater(a.t - b.t)(0,0) : g0_lesser(a.t - b.t)(0,0));
 }
};
