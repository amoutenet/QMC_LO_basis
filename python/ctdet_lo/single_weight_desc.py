# Generated automatically using the command :
# c++2py ../../c++/single_weight.hpp -C pytriqs --cxxflags="-std=c++17 -DHAS_OPTIONAL_HEADER "
from cpp2py.wrap_generator import *

# The module
module = module_(full_name = "single_weight", doc = "", app_name = "single_weight")

# Imports
import pytriqs.gf

# Add here all includes
module.add_include("single_weight.hpp")

# Add here anything to add in the C++ code at the start, e.g. namespace using
module.add_preamble("""
#include <cpp2py/converters/complex.hpp>
#include <triqs/cpp2py_converters/arrays.hpp>
#include <triqs/cpp2py_converters/gf.hpp>

""")


module.add_function ("double single_weight (array<double,1> times, array<int,1> i_taus, array<int,1> ls, double U, gf<triqs::gfs::retime> g0_lesser, gf<triqs::gfs::retime> g0_greater, triqs::gfs::dcomplex diag, triqs::gfs::dcomplex alpha, double t_max)", doc = """""")



module.generate_code()
