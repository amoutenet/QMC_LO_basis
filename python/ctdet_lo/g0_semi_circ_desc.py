# Generated automatically using the command :
# c++2py ../../c++/g0_semi_circ.hpp -C pytriqs --cxxflags="-std=c++17 -DHAS_OPTIONAL_HEADER "
from cpp2py.wrap_generator import *

# The module
module = module_(full_name = "g0_semi_circ", doc = "", app_name = "g0_semi_circ")

# Imports
import pytriqs.gf

# Add here all includes
module.add_include("g0_semi_circ.hpp")

# Add here anything to add in the C++ code at the start, e.g. namespace using
module.add_preamble("""
#include <cpp2py/converters/pair.hpp>
#include <triqs/cpp2py_converters/gf.hpp>

""")


module.add_function ("std::pair<gf_view<refreq>,gf_view<refreq>> make_g0_semi_circular (double beta, double Gamma, double tmax_gf0, int Nt_gf0, double epsilon_d)", doc = """""")



module.generate_code()