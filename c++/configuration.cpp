#include "./configuration.hpp"
#include "./parameters.hpp"

//#define CHECK_GRAY_CODE_INTEGRITY
//#define REGENERATE_MATRIX_BEFORE_EACH_GRAY_CODE

/// Gray code determinant rotation. Returns the sum of prod of det for all keldysh configurations.
dcomplex configuration::recompute_sum_keldysh_indices(std::vector<double> V) {

 if (order == 0) return imag(diag);

 dcomplex res{0,0};
 int sign = -1;                    // Starting with a flip, so -1 -> 1, which is needed in the first iteration
 auto two_to_k = uint64_t(1) << order; // shifts the bits k time to the left

 std::vector<keldysh_pt> vertices(V.size());
 for (int i = 0; i < V.size(); i++) vertices[i] = keldysh_pt{V[i],0};

 for (uint64_t n = 0; n < two_to_k; ++n) {
  int nlc = (n < two_to_k - 1 ? ffs(~n) : order) - 1;
  auto p = flip_index(vertices[nlc]);
  vertices[nlc] = p;

  res += sign * compute_det(0, vertices) * compute_det(1, vertices);
  sign = -sign;

 }

 dcomplex i_n[4] = {{1, 0}, {0, 1}, {-1, 0}, {0, -1}}; // powers of i
 res = res * -i_n[(order + 1) % 4];                        // * - i^(k+1)

 return real(res);
}

dcomplex configuration::compute_det(int sigma, std::vector<lo_vertex> V) {
  int dim = V.size() + (1 - sigma);
     
  if (sigma == 0) A(V.size(), V.size()) = g0_lo[sigma](vmax, vmax);
  for (int i =0; i < V.size(); i++) {
    if (sigma == 0) {
      A(V.size(), i) = g0_lo[sigma](vmax, V[i]);
      A(i, V.size()) = g0_lo[sigma](V[i], vmax);
    }
      
    for (int j =0; j < V.size(); j++) A(i,j) = g0_lo[sigma](V[i], V[j]);
  }
    
  return detclass.determinant(A, ri, dim);
}

dcomplex configuration::compute_det(int sigma, std::vector<keldysh_pt> V) {
  int dim = V.size() + (1 - sigma);
  
  if (sigma == 0) A(V.size(), V.size()) = g0_pm[sigma]({tmax,0}, {tmax,1});
    for (int i =0; i < V.size(); i++) {
      if (sigma == 0) {
	A(V.size(), i) = g0_pm[sigma]({tmax,0}, V[i]);
	A(i, V.size()) = g0_pm[sigma](V[i], {tmax,1});
      }

      for (int j =0; j < V.size(); j++) A(i,j) = g0_pm[sigma](V[i], V[j]);
    }
   
  return detclass.determinant(A, ri, dim);
}


std::pair<int, bool> configuration::get_flip_idx(std::vector<lo_vertex> V) { 
  int last_size = 0;
  int last_end = 0;
  bool last_all_l_equal = true;
	       
  int i = 0;
  while (i < V.size()) {
    auto itau = V[i].itau;
    int size = 1;

    for (int j = i+1; (j < V.size()) and (V[j].itau == itau); j++) size += 1;
    int end = i + size - 1;

    bool all_l_equal = true;
    auto l = V[i].l;
    for (int j = i + 1; (j < i + size - 1) && all_l_equal; j++) all_l_equal = (V[j].l == l);

    if ((size > last_size) || ((size == last_size) && (!all_l_equal || last_all_l_equal))) {
      last_size = size;
      last_end = end;
      last_all_l_equal = all_l_equal;
    }

    i += size;
  }

  return {last_end, last_all_l_equal};
}
