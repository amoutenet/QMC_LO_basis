#include <iostream>
#include <triqs/gfs.hpp>
#include "solver_core.hpp"
#include "g0_semi_circ.hpp"
#include "parameters.hpp"


using namespace triqs::gfs;
using namespace triqs::arrays;

int main(int argc, char* argv[]) {
 mpi::communicator world;
 mpi::environment env(argc, argv);

 solver_core s = solver_core();

 std::pair<gf<refreq>, gf<refreq>> g0 = make_g0_semi_circular(100, 0.2 * 0.2, 400, 20000, -0.36);

 solve_parameters_t params;
 params.measure = "n_lo";

 params.lambda = 0.0288;
 params.U = 1.2;
 params.tmax = 400;
 params.g0_R = std::get<0>(g0);
 params.g0_K = std::get<1>(g0);

 params.perturbation_order = 6;
 params.min_perturbation_order = 5;

 params.alpha = 0.3;
 params.cauchy_t = 380;
 params.cauchy_a = 37;

 params.n_cycles = 5000000;
 params.length_cycle = 14;
 params.n_warmup_cycles = 20000;

 s.solve(params);

 return 0;
}




