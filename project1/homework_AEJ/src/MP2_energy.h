/**
 * @file MP2_energy.h
 * @brief Declaration of helper functions for integral retrieval.
 *
 * This header defines the interface for accessing two-electron integrals (ERIs)
 * stored in sparse formats, handling the necessary symmetry index permutations.
 *
 * @date 2025
 */

#include <stdint.h>

/**
 * @brief Retrieves a specific Two-Electron Integral (ERI) from sparse arrays.
 *
 * This function searches for the integral \f$ (pq|rs) \f$ within the provided
 * sparse lists (index and value).
 *
 * **Algorithm details:**
 * Since integrals are often stored exploiting 8-fold permutational symmetry
 * (e.g., \f$ (pq|rs) = (rs|pq) = (qp|rs) \dots \f$), this function:
 * 1. Canonicalizes the input indices `p, q, r, s`.
 * 2. Iterates through the `index` array to find the match.
 * 3. Returns the corresponding `value` if found, or 0.0 otherwise.
 *
 * @param p Index of the first orbital.
 * @param q Index of the second orbital.
 * @param r Index of the third orbital.
 * @param s Index of the fourth orbital.
 * @param nint Total number of non-zero integrals stored (buffer size).
 * @param index Pointer to the flat array containing index quadruplets (size: \c 4*nint).
 * @param value Pointer to the array containing integral values (size: \c nint).
 * @return double The value of the integral \f$ (pq|rs) \f$.
 *
 * @code
 * double get_eri(int p, int q, int r, int s, int64_t nint, const int32_t* index, const double* value);
 * @endcode
 */
double get_eri(int p, int q, int r, int s, int64_t nint, const int32_t* index, const double* value);
