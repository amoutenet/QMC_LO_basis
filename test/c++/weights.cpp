#include <iostream>
#include <triqs/test_tools/gfs.hpp>
#include "weight_analysis.hpp"
#include "g0_semi_circ.hpp"

using namespace triqs::gfs;

TEST(Ctdet_LO, weights) {
 mpi::communicator world;
 std::pair<gf<refreq>, gf<refreq>> g0 = make_g0_semi_circular(10, 0.5 * 0.5, 100, 1000, 0.2);
  
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
 
 array<double, 1> times{380, 385, 390};
 auto weights = weight_analysis(times, 1.2, g0_lesser_t, g0_greater_t, diag, 0.3, 400).second;

 // Save the results
 if (world.rank() == 0) {
  triqs::h5::file G_file("weights.out.h5", 'w');
  h5_write(G_file, "weights", weights);
 }

 array<double, 1> w;

 if (world.rank() == 0) {
  triqs::h5::file G_file("weights.ref.h5", 'r');
  h5_read(G_file, "weights", w);
  EXPECT_ARRAY_NEAR(weights, w);
 }
}

MAKE_MAIN
