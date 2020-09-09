/*******************************************************************************
 *
 * Computes the determinant of the n x n matrix A. It is obtained from
 * a product of the diagonal elements of U in the factorization of
 * the matrix A into A = LU where L is lower triangular and U is upper
 * triangular. The matrix is overwritten by the LU decomposition.
 * The diagonal elements of L (which are unity) are not stored.
 *
 * Copyright (C) 2018 by M. Ferrero
 *
 ******************************************************************************/
#pragma once

#include <triqs/arrays.hpp>
#include <vector>

struct michel_determinant {

  michel_determinant() {}

std::complex<double> determinant(triqs::arrays::array_view<std::complex<double>,2> A, std::vector<int> & ri, int n) {

  std::complex<double> det = 1.0;

  // Initialize the pointer vector
  for (int i = 0; i < n; i++) ri[i] = i;

  // LU factorization
  for (int p = 1; p < n; p++) {

    // Find pivot element in column [p-1]
    for (int i = p; i < n; i++) {
      if (abs(A(ri[i], p-1)) > abs(A(ri[p-1], p-1))) {
        // Switch the index for the p-1 pivot row if necessary
        int t = ri[p-1];
        ri[p-1] = ri[i];
        ri[i] = t;
        det = -det;
      }
    }

    // Multiply the diagonal elements
    det *= A(ri[p-1], p-1);

    // Form multiplier
    auto coef = 1 / (A(ri[p-1], p-1));
    for (int i = p; i < n; i++) {
      A(ri[i], p-1) *= coef;
      // Eliminate [p-1]
      for (int j = p; j < n; j++)
        A(ri[i], j) -= A(ri[i], p-1) * A(ri[p-1], j);
    }

  }

  det *= A(ri[n-1], n-1);
  return det;

}
};
