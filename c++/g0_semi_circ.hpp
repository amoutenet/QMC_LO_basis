#pragma once
#include <triqs/gfs.hpp>
#include <triqs/arrays.hpp>

using namespace triqs::gfs;

// Semi-circular band
std::pair<gf<refreq>, gf<refreq>> make_g0_semi_circular(double beta, double Gamma, double tmax_gf0, int Nt_gf0,
                                                                  double epsilon_d) {

 // Construction of the empty GF's, with the correct number of points
 auto g0_R_t = gf<retime>{{-tmax_gf0, tmax_gf0, 2 * Nt_gf0}, {2, 2}}; // even number of points needed -- see solver_core
 auto g0_R_w = make_gf_from_fourier(g0_R_t, true); // we want a point at 0 in frequency

 auto g0_K_t = g0_R_t;
 auto g0_K_w = g0_R_w;

 // Retarded self energy for the semi-circular bath
 auto sigma_bath = [](double omega) -> dcomplex {
  omega = omega / 2.;
  if (std::abs(omega) < 1) return dcomplex{omega, -std::sqrt(1 - omega * omega)};
  if (omega > 1) return omega - std::sqrt(omega * omega - 1);
  return omega + std::sqrt(omega * omega - 1);
 };

 // The non interacting dot GF's in frequency (2*2 matrix with Keldysh indices)
 auto G0_dd_w = [&](double omega) {
  dcomplex gr = 1 / (omega - epsilon_d - Gamma * sigma_bath(omega));
  dcomplex ga = conj(gr);
  dcomplex gk = std::tanh(beta * omega / 2.) * (gr - ga);
  return array<dcomplex, 2>{{gr, gk}, {dcomplex{0,0}, ga}};

  //dcomplex fac = 2_j * imag(gr);
  //dcomplex l = - nf(omega) * fac; // G+- = G lesser
  //dcomplex g = l + fac; // G-+ = G greater // to be replaced by KMS
  //return array<dcomplex, 2>{{gr + l, l}, {g, -conj(gr) + l}};
  // 0,0: G++ = GR  + G+-/ 0,1: G+- = G</ 1,0: G-+ = G>/ 1,1: G-- = G+- - GA = G+- - (GR)*
 };

 for (auto w : g0_R_w.mesh()) {
  auto g0_dd = G0_dd_w(w);
  g0_R_w[w](0, 0) = g0_dd(0, 0);
  g0_K_w[w](0, 0) = g0_dd(0, 1);
  // Never used for density
  g0_R_w[w](0, 1) = 0.0; //g0_dc(0, 1);
  g0_K_w[w](0, 1) = 0.0; //g0_dc(1, 0);
  // Set lower components to zero
  g0_R_w[w](1, 0) = 0.0;
  g0_R_w[w](1, 1) = 0.0;
  g0_K_w[w](1, 0) = 0.0;
  g0_K_w[w](1, 1) = 0.0;
 }

 // No singularity information needed.

 return {g0_R_w, g0_K_w};
}
