FILE(GLOB_RECURSE cpp_sources RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.cpp)
add_library(ctdet_lo_c ${cpp_sources})
target_link_libraries(ctdet_lo_c PUBLIC triqs)

install(TARGETS ctdet_lo_c DESTINATION lib)
