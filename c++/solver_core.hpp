#include <triqs/gfs.hpp>
#include <triqs/arrays.hpp>
#include <triqs/mc_tools/mc_generic.hpp>
#include <triqs/statistics/histograms.hpp>

#include "parameters.hpp"

// ------------ The main class of the solver -----------------------

using namespace triqs::gfs;

class solver_core {
 //mpi communicator
 mpi::communicator world;

 public:

 triqs::arrays::vector<double> d;
 triqs::arrays::vector<double> eta;
 triqs::arrays::vector<double> proba;
 triqs::arrays::vector<double> av_sign;

 std::optional<triqs::statistics::histogram> times_histogram; // Histogram of the visited times at the highest perturbation order
 std::optional<triqs::statistics::histogram> weights_histogram; // Histogram of the MC weights at the highest perturbation order
 std::optional<triqs::arrays::vector<double>> correlation;
 std::optional<triqs::arrays::matrix<double>> times;

 triqs::gfs::gf<triqs::gfs::retime> g0_lesser_t;
 triqs::gfs::gf<triqs::gfs::retime> g0_greater_t;

 triqs::gfs::gf<triqs::gfs::refreq> g0_lesser_w;
 triqs::gfs::gf<triqs::gfs::refreq> g0_greater_w;

 solver_core(): d(1), eta(1), proba(1), av_sign(1) {}

 CPP2PY_ARG_AS_DICT // Wrap the solver parameters as a ** call in python with the clang & c++2py tool
 void solve(solve_parameters_t const &params);
};
