/**
 * @file MP2_energy.c
 * @brief Implementation of helper functions for MP2 energy calculation.
 *
 * This file contains the logic to handle the retrieval of Two-Electron Integrals (ERIs)
 * from the sparse data structures provided by TREXIO.
 *
 * @author Asia Meli, Yiyi Yang and Eloá Abreu
 * @date 2025
**/

#include "MP2_energy.h"
/**
 * @brief Retrieve a two-electron integral (ERI) from TREXIO sparse storage using 8-fold permutational symmetry.
 *
 * TREXIO stores electron-repulsion integrals (ERIs) in a sparse format:
 * an array of indices (i,j,k,l) and a corresponding array of values v = <ij|kl>.
 * Since ERIs satisfy 8-fold permutational symmetry, the same numerical value can be
 * represented by any of the following equivalent index orderings:
 *
 * \f[
 * (ij|kl) = (ji|lk) = (il|kj) = (li|jk) = (kl|ij) = (lk|ji) = (kj|il) = (jk|li)
 * \f]
 *
 * This function performs a linear scan over the list of stored integrals and returns
 * the value if any of the 8 symmetric permutations matches the requested (p,q,r,s).
 * If no match is found in the current buffer, the function returns 0.0.
 *
 * @param p First orbital index.
 * @param q Second orbital index.
 * @param r Third orbital index.
 * @param s Fourth orbital index.
 * @param nint Number of integrals stored in the buffer (length of @p value, and @p index/4).
 * @param index Pointer to the flattened index array of size 4*nint:
 *              index[4*n + 0..3] = {i, j, k, l} for the n-th integral.
 * @param value Pointer to the integral values array of size nint.
 *
 * @return double The ERI value for the requested indices (including symmetry matches),
 *         or 0.0 if the integral is not present in the sparse list.
 *
 * @note Complexity is O(nint) per query because of the linear scan. This is acceptable
 *       for small systems or demonstration purposes, but can be expensive for large datasets.
 */
  

double get_eri(int p, int q, int r, int s, int64_t nint, const int32_t* index, const double* value) {

  for (int64_t n = 0; n < nint; n++) {
    int i = index[4*n + 0];
    int j = index[4*n + 1];
    int k = index[4*n + 2];
    int l = index[4*n + 3];
    double v = value[n];

    if (
      (i==p && j==q && k==r && l==s) ||  // <ij|kl>
      (i==p && l==q && k==r && j==s) ||  // <il|kj>
      (k==p && l==q && i==r && j==s) ||  // <kl|ij>
      (k==p && j==q && i==r && l==s) ||  // <kj|il>
      (j==p && i==q && l==r && k==s) ||  // <ji|lk>
      (l==p && i==q && j==r && k==s) ||  // <li|jk>
      (l==p && k==q && j==r && i==s) ||  // <lk|ji>
      (j==p && k==q && l==r && i==s)     // <jk|li>
    ) {
      return v;
    }
  }
// If not found in the list, the integral is zero (sparse format) 
  return 0.0;
}
