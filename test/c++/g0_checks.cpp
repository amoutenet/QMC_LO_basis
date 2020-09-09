#include <iostream>
#include <triqs/test_tools/gfs.hpp>
#include "g0_semi_circ.hpp"

using namespace triqs::gfs;

TEST(Ctdet_LO, g0_checks) {
 mpi::communicator world;
 std::pair<gf<refreq>, gf<refreq>> g0 = make_g0_semi_circular(10, 0.5 * 0.5, 100, 1000, 0.2);

 // Save the results
 if (world.rank() == 0) {
  triqs::h5::file G_file("g0_checks.out.h5", 'w');
  h5_write(G_file, "G_R", std::get<0>(g0));
  h5_write(G_file, "G_K", std::get<1>(g0));
 }

 gf<refreq> g_r;
 gf<refreq> g_k;

 if (world.rank() == 0) {
  triqs::h5::file G_file("g0_checks.ref.h5", 'r');
  h5_read(G_file, "G_R", g_r);
  h5_read(G_file, "G_K", g_k);
  EXPECT_GF_NEAR(std::get<0>(g0), g_r);
  EXPECT_GF_NEAR(std::get<1>(g0), g_k);
 }
}

MAKE_MAIN
