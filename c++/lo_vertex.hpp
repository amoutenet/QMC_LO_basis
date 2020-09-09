#pragma once
#include <triqs/mc_tools/random_generator.hpp>

// ------------------- Vertex in the LO basis -------------------------
/// A LO vertex. It is described by its time, its interaction matrix 0.5 * (1*tau + tau*1), and an index l for the identity.
struct lo_vertex{
  double t; // time, in [0, t_max].
  int itau; // itau = sigma: spin sigma has vertex tau
  int l; // l = 0 or 1, LO index associated to the identity

};

inline  bool operator<(lo_vertex const & v1, lo_vertex const & v2) {
    return v1.t < v2.t;
  }

inline lo_vertex get_random_lo_vertex(double t_max, triqs::mc_tools::random_generator* RNG, double cauchy_t, double cauchy_a) {
 /// The following lines generate a u distributed according a Cauchy law - see notes
 double t = (*RNG)(1.);
 auto x = (1 - t)* std::atan(- cauchy_t / cauchy_a) + t * std::atan((t_max - cauchy_t) / cauchy_a);
 auto u = cauchy_t + cauchy_a * std::tan(x);

 return lo_vertex{u, (*RNG)(2), (*RNG)(2)};
}
