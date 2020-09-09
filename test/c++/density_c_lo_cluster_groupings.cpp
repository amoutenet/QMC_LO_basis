#include <iostream>
#include <triqs/gfs.hpp>
#include "solver_core.hpp"
#include "g0_semi_circ.hpp"
#include "parameters.hpp"
#include <triqs/test_tools/arrays.hpp>


using namespace triqs::gfs;
using namespace triqs::arrays;

TEST(Ctdet_LO, density_c_lo_cluster_groupings) {
 mpi::communicator world;

 solver_core s = solver_core();

 std::pair<gf<refreq>, gf<refreq>> g0 = make_g0_semi_circular(100, 0.2 * 0.2, 400, 20000, -0.36);

 solve_parameters_t params;
 params.measure = "n_lo_cluster_groupings";
 params.U = 1.2;
 params.tmax = 400;
 params.g0_R = std::get<0>(g0);
 params.g0_K = std::get<1>(g0);
 params.perturbation_order = 3;
 params.min_perturbation_order = 3;

 params.alpha = 0.3;
 params.cauchy_t = 379;
 params.cauchy_a = 23;

 params.n_cycles = 30;
 params.length_cycle = 10;
 params.n_warmup_cycles = 5;

 params.lambda = 0.0387;

 s.solve(params);

 if (world.rank() == 0) {
  triqs::h5::file G_file("density_c_lo_cluster_groupings.out.h5", 'w');
  h5_write(G_file, "d", s.d);
  h5_write(G_file, "eta", s.eta);
  h5_write(G_file, "proba", s.proba);
 }

 array<double, 1> d;
 array<double, 1> eta;
 array<double, 1> proba;
 array<double, 1> av_sign;

 if (world.rank() == 0) {
  triqs::h5::file G_file("density_c_lo_cluster_groupings.ref.h5", 'r');
  h5_read(G_file, "d", d);
  h5_read(G_file, "eta", eta);
  h5_read(G_file, "proba", proba);

  const double pre = 1e-12;
  EXPECT_ARRAY_NEAR(d, s.d, pre);
  EXPECT_ARRAY_NEAR(eta, s.eta, pre);
  EXPECT_ARRAY_NEAR(proba, s.proba, pre);
 }

}

MAKE_MAIN
