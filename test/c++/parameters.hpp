#pragma once
#include <triqs/gfs.hpp>

using namespace triqs::utility;

// All the arguments of the solve function
struct solve_parameters_t {
 /// lambda parameter - used to adjust the time spent in each order
 double lambda;
 /// physical U - needed for the alpha shift
 double U;

 /// tmax
 double tmax;

 /// g0_R and g0_K
 triqs::gfs::gf<triqs::gfs::refreq> g0_R;
 triqs::gfs::gf<triqs::gfs::refreq> g0_K;

 /// alpha term
 double alpha = 0;

 /// Measure: density n_mix (mix algorithm), n_lo (pure LO code), n_mix_cancel_sym, n_lo_cancel_sym, MC times
 std::string measure = "n_lo";
 
 /// In case we measure MC times, the size of the times array
 long n_times = 100000;

 /// Adjusting the choice of times with a Cauchy law, because of clusterisation
 double cauchy_t;
 double cauchy_a;

 /// Do we want time and weight histograms?
 bool measure_histos = false;
 /// Do we want to measure decoorelation?
 bool measure_correlation = false;

 // --- Haar scattering parameters
 std::vector<triqs::arrays::matrix<int>> perm_vecs; //= std::vector<triqs::arrays::matrix<int>>{{0}};
 //std::vector<std::vector<int>> perm_vec = std::vector<std::vector<int>>{{0}};
 //std::vector<std::vector<int>> min_perm_vec = std::vector<std::vector<int>>{{0}};

 // ----   QMC parameters

 /// Perturbation order in U we mainly want to sample
 int perturbation_order;
 /// Minimal perturbation order the algm can reach (by default zero)
 int min_perturbation_order = 0;

 /// Number of QMC cycles
 long n_cycles;

 /// Length of a single QMC cycle
 int length_cycle;

 /// Number of cycles for thermalization
 int n_warmup_cycles;

 /// Seed for random number generator
 int random_seed = 34788 + 928374 * mpi::communicator().rank();

 /// Name of random number generator
 std::string random_name = "";

 /// Maximum runtime in seconds, use -1 to set infinite
 int max_MC_time = -1;

 /// Verbosity level
 int verbosity = ((mpi::communicator().rank() == 0) ? 3 : 0); // silence the slave nodes
};
