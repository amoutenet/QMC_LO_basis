# Build the python module
add_cpp2py_module(solver_core)
add_cpp2py_module(g0_semi_circ)
add_cpp2py_module(weight_analysis)
add_cpp2py_module(haar_analysis)
add_cpp2py_module(single_weight)

target_link_libraries(solver_core ctdet_lo_c)
target_link_libraries(g0_semi_circ ctdet_lo_c)
target_link_libraries(weight_analysis ctdet_lo_c)
target_link_libraries(haar_analysis ctdet_lo_c)
target_link_libraries(single_weight ctdet_lo_c)

# We need to include the convertes.hxx files
target_include_directories(solver_core PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

# Add version and __init__ file
configure_file(version.py.in version.py)
configure_file(__init__.py __init__.py COPYONLY)

set(PYTHON_SOURCES
 ${CMAKE_CURRENT_SOURCE_DIR}/__init__.py
 ${CMAKE_CURRENT_BINARY_DIR}/version.py
 )

# Install python module to proper location
set(PYTHON_LIB_DEST ${CPP2PY_PYTHON_LIB_DEST_ROOT}/triqs_ctdet_lo)
install(TARGETS solver_core DESTINATION ${PYTHON_LIB_DEST})
install(TARGETS g0_semi_circ DESTINATION ${PYTHON_LIB_DEST})
install(TARGETS weight_analysis DESTINATION ${PYTHON_LIB_DEST})
install(TARGETS haar_analysis DESTINATION ${PYTHON_LIB_DEST})
install(TARGETS single_weight DESTINATION ${PYTHON_LIB_DEST})
install(FILES ${PYTHON_SOURCES} DESTINATION ${PYTHON_LIB_DEST})

