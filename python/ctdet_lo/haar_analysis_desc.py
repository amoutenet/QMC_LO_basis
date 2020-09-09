# Generated automatically using the command :
# c++2py ../../c++/haar_analysis.hpp -C pytriqs --cxxflags="-std=c++17 -DHAS_OPTIONAL_HEADER "
from cpp2py.wrap_generator import *

# The module
module = module_(full_name = "haar_analysis", doc = "", app_name = "haar_analysis")

# Imports
import pytriqs.gf

# Add here all includes
module.add_include("haar_analysis.hpp")

# Add here anything to add in the C++ code at the start, e.g. namespace using
module.add_preamble("""
#include <cpp2py/converters/complex.hpp>
#include <triqs/cpp2py_converters/arrays.hpp>
#include <triqs/cpp2py_converters/gf.hpp>

""")


module.add_function ("std::pair<array<int, 1>, array<double,1>> haar_analysis (array<double,1> times, double U, gf<triqs::gfs::retime> g0_lesser, gf<triqs::gfs::retime> g0_greater, triqs::gfs::dcomplex diag, triqs::gfs::dcomplex alpha, double t_max)", doc = """""")



module.generate_code()
