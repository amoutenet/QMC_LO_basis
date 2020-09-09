#include <iostream>
#include <triqs/gfs.hpp>
#include <triqs/mc_tools/random_generator.hpp>
#include <configuration.hpp>
#include "solver_core.hpp"
#include "g0_semi_circ.hpp"
#include "lo_vertex.hpp"
#include <triqs/test_tools/gfs.hpp>

using namespace triqs::gfs;

TEST(Ctdet_LO, moves) {
 mpi::communicator world;

 solver_core s = solver_core();

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

 configuration config = configuration(g0_lesser_t, g0_greater_t, 0, 50, diag, 1);

 triqs::mc_tools::random_generator RNG("mt19937", 23432);

 std::vector<lo_vertex> V(config.order + 1);
 
 for (int i =0; i < V.size(); i++) V[i] = get_random_lo_vertex(50, &RNG, 48.5, 2.3);

 std::sort(V.begin(), V.end());

 std::complex<double> det = 1.0;
 // Construt det up and det down
 for (int sigma = 0; sigma < 2; sigma++) det *= config.compute_det(sigma, V);

 // Save the results
 if (world.rank() == 0) {
  triqs::h5::file G_file("moves.out.h5", 'w');
  h5_write(G_file, "det", det);
 }

 dcomplex d;

 if (world.rank() == 0) {
  triqs::h5::file G_file("moves.ref.h5", 'r');
  h5_read(G_file, "det", d);

  const double pre = 1e-12;
  EXPECT_COMPLEX_NEAR(det, d, pre);
 }
}

MAKE_MAIN
