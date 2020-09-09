# Generated automatically using the command :
# c++2py ../../c++/solver_core.hpp -p -msolver_core -o solver_core --moduledoc="The ctint solver" -C pytriqs --cxxflags="-std=c++17 -DHAS_OPTIONAL_HEADER -stdlib=libc++" --cxxflags="-std=c++17"
from cpp2py.wrap_generator import *

# The module
module = module_(full_name = "solver_core", doc = "The ctint solver", app_name = "solver_core")

# Imports
import pytriqs.gf
import pytriqs.statistics.histograms

# Add here all includes
module.add_include("solver_core.hpp")

# Add here anything to add in the C++ code at the start, e.g. namespace using
module.add_preamble("""
#include <cpp2py/converters/optional.hpp>
#include <cpp2py/converters/string.hpp>
#include <cpp2py/converters/vector.hpp>
#include <triqs/cpp2py_converters/arrays.hpp>
#include <triqs/cpp2py_converters/gf.hpp>

""")


# The class solver_core
c = class_(
        py_type = "SolverCore",  # name of the python class
        c_type = "solver_core",   # name of the C++ class
        doc = """""",   # doc of the C++ class
        hdf5 = False,
)

c.add_member(c_name = "d",
             c_type = "triqs::arrays::vector<double>",
             read_only= False,
             doc = """""")

c.add_member(c_name = "eta",
             c_type = "triqs::arrays::vector<double>",
             read_only= False,
             doc = """""")

c.add_member(c_name = "proba",
             c_type = "triqs::arrays::vector<double>",
             read_only= False,
             doc = """""")

c.add_member(c_name = "av_sign",
             c_type = "triqs::arrays::vector<double>",
             read_only= False,
             doc = """""")

c.add_member(c_name = "times_histogram",
             c_type = "std::optional<triqs::statistics::histogram>",
             read_only= False,
             doc = """""")

c.add_member(c_name = "weights_histogram",
             c_type = "std::optional<triqs::statistics::histogram>",
             read_only= False,
             doc = """""")

c.add_member(c_name = "correlation",
             c_type = "std::optional<triqs::arrays::vector<double> >",
             read_only= False,
             doc = """""")

c.add_member(c_name = "times",
             c_type = "std::optional<triqs::arrays::matrix<double> >",
             read_only= False,
             doc = """""")

c.add_member(c_name = "g0_lesser_t",
             c_type = "triqs::gfs::gf<triqs::gfs::retime>",
             read_only= False,
             doc = """""")

c.add_member(c_name = "g0_greater_t",
             c_type = "triqs::gfs::gf<triqs::gfs::retime>",
             read_only= False,
             doc = """""")

c.add_member(c_name = "g0_lesser_w",
             c_type = "triqs::gfs::gf<triqs::gfs::refreq>",
             read_only= False,
             doc = """""")

c.add_member(c_name = "g0_greater_w",
             c_type = "triqs::gfs::gf<triqs::gfs::refreq>",
             read_only= False,
             doc = """""")

c.add_constructor("""()""", doc = """""")

c.add_method("""void solve (**solve_parameters_t)""",
             doc = """
+------------------------+------------------------------------+------------------------------------------------+-----------------------------------------------------------------------------------------------------------+
| Parameter Name         | Type                               | Default                                        | Documentation                                                                                             |
+========================+====================================+================================================+===========================================================================================================+
| lambda                 | double                             |                                                | lambda parameter - used to adjust the time spent in each order                                            |
+------------------------+------------------------------------+------------------------------------------------+-----------------------------------------------------------------------------------------------------------+
| U                      | double                             |                                                | physical U - needed for the alpha shift                                                                   |
+------------------------+------------------------------------+------------------------------------------------+-----------------------------------------------------------------------------------------------------------+
| tmax                   | double                             |                                                | tmax                                                                                                      |
+------------------------+------------------------------------+------------------------------------------------+-----------------------------------------------------------------------------------------------------------+
| g0_R                   | triqs::gfs::gf<triqs::gfs::refreq> |                                                | g0_R and g0_K                                                                                             |
+------------------------+------------------------------------+------------------------------------------------+-----------------------------------------------------------------------------------------------------------+
| g0_K                   | triqs::gfs::gf<triqs::gfs::refreq> |                                                |                                                                                                           |
+------------------------+------------------------------------+------------------------------------------------+-----------------------------------------------------------------------------------------------------------+
| alpha                  | double                             | 0                                              | alpha term                                                                                                |
+------------------------+------------------------------------+------------------------------------------------+-----------------------------------------------------------------------------------------------------------+
| measure                | std::string                        | "n_mix"                                        | Measure: density n_mix (mix algorithm), n_lo (pure LO code), n_mix_cancel_sym, n_lo_cancel_sym, MC times  |
+------------------------+------------------------------------+------------------------------------------------+-----------------------------------------------------------------------------------------------------------+
| n_times                | long                               | 100000                                         | In case we measure MC times, the size of the times array                                                  |
+------------------------+------------------------------------+------------------------------------------------+-----------------------------------------------------------------------------------------------------------+
| cauchy_t               | double                             |                                                | Adjusting the choice of times with a Cauchy law, because of clusterisation                                |
+------------------------+------------------------------------+------------------------------------------------+-----------------------------------------------------------------------------------------------------------+
| cauchy_a               | double                             |                                                |                                                                                                           |
+------------------------+------------------------------------+------------------------------------------------+-----------------------------------------------------------------------------------------------------------+
| measure_histos         | bool                               | false                                          | Do we want time and weight histograms?                                                                    |
+------------------------+------------------------------------+------------------------------------------------+-----------------------------------------------------------------------------------------------------------+
| measure_correlation    | bool                               | false                                          | Do we want to measure decoorelation?                                                                      |
+------------------------+------------------------------------+------------------------------------------------+-----------------------------------------------------------------------------------------------------------+
| perm_vecs              | std::vector<triqs::arrays::matrix<int>>      |                                      |                                                                                                           |
+------------------------+------------------------------------+------------------------------------------------+-----------------------------------------------------------------------------------------------------------+
| perturbation_order     | int                                |                                                | Perturbation order in U we mainly want to sample                                                          |
+------------------------+------------------------------------+------------------------------------------------+-----------------------------------------------------------------------------------------------------------+
| min_perturbation_order | int                                | 0                                              | Minimal perturbation order the algm can reach (by default zero)                                           |
+------------------------+------------------------------------+------------------------------------------------+-----------------------------------------------------------------------------------------------------------+
| n_cycles               | long                               |                                                | Number of QMC cycles                                                                                      |
+------------------------+------------------------------------+------------------------------------------------+-----------------------------------------------------------------------------------------------------------+
| length_cycle           | int                                |                                                | Length of a single QMC cycle                                                                              |
+------------------------+------------------------------------+------------------------------------------------+-----------------------------------------------------------------------------------------------------------+
| n_warmup_cycles        | int                                |                                                | Number of cycles for thermalization                                                                       |
+------------------------+------------------------------------+------------------------------------------------+-----------------------------------------------------------------------------------------------------------+
| random_seed            | int                                | 34788+928374*mpi::communicator().rank() | Seed for random number generator                                                                          |
+------------------------+------------------------------------+------------------------------------------------+-----------------------------------------------------------------------------------------------------------+
| random_name            | std::string                        | ""                                             | Name of random number generator                                                                           |
+------------------------+------------------------------------+------------------------------------------------+-----------------------------------------------------------------------------------------------------------+
| max_MC_time            | int                                | -1                                             | Maximum runtime in seconds, use -1 to set infinite                                                        |
+------------------------+------------------------------------+------------------------------------------------+-----------------------------------------------------------------------------------------------------------+
| verbosity              | int                                | ((mpi::communicator().rank()==0)?3:0)   | Verbosity level                                                                                           |
+------------------------+------------------------------------+------------------------------------------------+-----------------------------------------------------------------------------------------------------------+
""")

module.add_class(c)


# Converter for solve_parameters_t
c = converter_(
        c_type = "solve_parameters_t",
        doc = """""",
)
c.add_member(c_name = "lambda",
             c_type = "double",
             initializer = """  """,
             doc = """lambda parameter - used to adjust the time spent in each order""")

c.add_member(c_name = "U",
             c_type = "double",
             initializer = """  """,
             doc = """physical U - needed for the alpha shift""")

c.add_member(c_name = "tmax",
             c_type = "double",
             initializer = """  """,
             doc = """tmax""")

c.add_member(c_name = "g0_R",
             c_type = "triqs::gfs::gf<triqs::gfs::refreq>",
             initializer = """  """,
             doc = """g0_R and g0_K""")

c.add_member(c_name = "g0_K",
             c_type = "triqs::gfs::gf<triqs::gfs::refreq>",
             initializer = """  """,
             doc = """""")

c.add_member(c_name = "alpha",
             c_type = "double",
             initializer = """ 0 """,
             doc = """alpha term""")

c.add_member(c_name = "measure",
             c_type = "std::string",
             initializer = """ "n_mix" """,
             doc = """Measure: density n_mix (mix algorithm), n_lo (pure LO code), n_mix_cancel_sym, n_lo_cancel_sym, MC times""")

c.add_member(c_name = "n_times",
             c_type = "long",
             initializer = """ 100000 """,
             doc = """In case we measure MC times, the size of the times array""")

c.add_member(c_name = "cauchy_t",
             c_type = "double",
             initializer = """  """,
             doc = """Adjusting the choice of times with a Cauchy law, because of clusterisation""")

c.add_member(c_name = "cauchy_a",
             c_type = "double",
             initializer = """  """,
             doc = """""")

c.add_member(c_name = "measure_histos",
             c_type = "bool",
             initializer = """ false """,
             doc = """Do we want time and weight histograms?""")

c.add_member(c_name = "measure_correlation",
             c_type = "bool",
             initializer = """ false """,
             doc = """Do we want to measure decoorelation?""")

c.add_member(c_name = "perm_vecs",
        c_type = "std::vector<triqs::arrays::matrix<int>>",
        initializer = """  """,
             doc = """ """)

c.add_member(c_name = "perturbation_order",
             c_type = "int",
             initializer = """  """,
             doc = """Perturbation order in U we mainly want to sample""")

c.add_member(c_name = "min_perturbation_order",
             c_type = "int",
             initializer = """ 0 """,
             doc = """Minimal perturbation order the algm can reach (by default zero)""")

c.add_member(c_name = "n_cycles",
             c_type = "long",
             initializer = """  """,
             doc = """Number of QMC cycles""")

c.add_member(c_name = "length_cycle",
             c_type = "int",
             initializer = """  """,
             doc = """Length of a single QMC cycle""")

c.add_member(c_name = "n_warmup_cycles",
             c_type = "int",
             initializer = """  """,
             doc = """Number of cycles for thermalization""")

c.add_member(c_name = "random_seed",
             c_type = "int",
             initializer = """ 34788+928374*mpi::communicator().rank() """,
             doc = """Seed for random number generator""")

c.add_member(c_name = "random_name",
             c_type = "std::string",
             initializer = """ "" """,
             doc = """Name of random number generator""")

c.add_member(c_name = "max_MC_time",
             c_type = "int",
             initializer = """ -1 """,
             doc = """Maximum runtime in seconds, use -1 to set infinite""")

c.add_member(c_name = "verbosity",
             c_type = "int",
             initializer = """ ((mpi::communicator().rank()==0)?3:0) """,
             doc = """Verbosity level""")

module.add_converter(c)


module.generate_code()
