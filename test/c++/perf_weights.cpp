#include <iostream>
#include <triqs/gfs.hpp>
#include "weight_analysis.hpp"
#include "g0_semi_circ.hpp"
#include "parameters.hpp"

using namespace triqs::gfs;
using namespace triqs::arrays;

int main(int argc, char* argv[]) {
 mpi::communicator world;
 mpi::environment env(argc, argv);

 std::pair<gf<refreq>, gf<refreq>> g0 = make_g0_semi_circular(10, 0.5 * 0.5, 100, 10000, 0.2);

 auto g0_lesser_w = std::get<0>(g0);
 auto g0_greater_w = g0_lesser_w;
 
 for (auto w : g0_lesser_w.mesh()) {
   dcomplex fac = 2_j * imag(std::get<0>(g0)[w](0, 0));
   g0_lesser_w[w](0, 0) = (std::get<1>(g0)[w](0,0) - fac) / 2;
   g0_greater_w[w](0, 0) = (std::get<1>(g0)[w](0,0) + fac) / 2;
 }

 auto g0_lesser_t = make_gf_from_fourier(g0_lesser_w);
 auto g0_greater_t = make_gf_from_fourier(g0_greater_w);

 auto diag = make_gf_from_fourier(g0_lesser_w, true)(0)(0,0);

 array<double, 1> times{360, 375, 379, 380, 388, 390};
 auto weights = weight_analysis(times, 1.2, g0_lesser_t, g0_greater_t, diag, 0.3, 400);
 
 return 0;
}




