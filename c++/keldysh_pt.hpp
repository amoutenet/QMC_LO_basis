#pragma once

#pragma once
#include <triqs/mc_tools/random_generator.hpp>

// --------------   Point on the Keldysh contour   --------------------------------------

/// A point in time on the double contour, with an additionnal index x_index_t
struct keldysh_pt {
 double t; // time, in [0, t_max].
 int k;  // Keldysh index : 0 (upper contour), or 1 (lower contour)
};

/// flip the Keldysh index
inline keldysh_pt flip_index(keldysh_pt const &t) { return {t.t, 1 - t.k}; }
