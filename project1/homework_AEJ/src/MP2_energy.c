/**
 * @file MP2_energy.c
 * @brief Implementation of helper functions for MP2 energy calculation.
 *
 * This file contains the logic to handle the retrieval of Two-Electron Integrals (ERIs)
 * from the sparse data structures provided by TREXIO.
 *
 * @author Asia Meli, Yiyi Yang and Eloá Abreu
 * @date 2025
 */

#include "MP2_energy.h"
#include <stdlib.h>
#include <stdio.h>

/**
 * @brief Retrieves a specific Two-Electron Integral (ERI) considering symmetry.
 *
 * This function handles the 8-fold permutational symmetry of real ERIs:
 * \f[ (pq|rs) = (qp|rs) = (pq|sr) = (qp|sr) = (rs|pq) = (sr|pq) = (rs|qp) = (sr|qp) \f]
 *
 * The algorithm works as follows:
 *
 * **Step 1: Define Permutations**
 * Since we don't know which canonical version is stored in the file, we prepare
 * all 8 possible combinations of the requested indices (p,q,r,s).
 * @code
 * int permutations[8][4] = {
 * {p, q, r, s}, {q, p, r, s},
 * {p, q, s, r}, {q, p, s, r},
 * {r, s, p, q}, {s, r, p, q},
 * {r, s, q, p}, {s, r, q, p}
 * };
 * @endcode
 *
 * **Step 2: Linear Search**
 * We iterate through the sparse list of integrals. For each stored integral,
 * we check if its indices match ANY of our 8 permutations.
 * @code
 * for (int64_t k = 0; k < nint; k++) {
 * // ... check if stored indices match any permutation ...
 * if (match_found) return value[k];
 * }
 * @endcode
 *
 * @param p Index of the first orbital.
 * @param q Index of the second orbital.
 * @param r Index of the third orbital.
 * @param s Index of the fourth orbital.
 * @param nint Total number of integrals in the buffer.
 * @param index Pointer to the flattened index array (size 4*nint).
 * @param value Pointer to the value array (size nint).
 * @return double The value of the integral if found, 0.0 otherwise.
 */
double get_eri(int p, int q, int r, int s, int64_t nint, const int32_t* index, const double* value) {

    // Array of the 8 possible permutations of indices (p,q,r,s)
    int permutations[8][4] = {
        {p, q, r, s},
        {q, p, r, s},
        {p, q, s, r},
        {q, p, s, r},
        {r, s, p, q},
        {s, r, p, q},
        {r, s, q, p},
        {s, r, q, p}
    };

    // Iterate through the sparse list of integrals
    for (int64_t k = 0; k < nint; k++) {

        // Fetch the indices stored at position k
        int i_stored = index[4 * k + 0];
        int j_stored = index[4 * k + 1];
        int k_stored = index[4 * k + 2];
        int l_stored = index[4 * k + 3];

        // Check against all 8 permutations
        for (int m = 0; m < 8; m++) {
            if (permutations[m][0] == i_stored &&
                permutations[m][1] == j_stored &&
                permutations[m][2] == k_stored &&
                permutations[m][3] == l_stored) {

                return value[k]; // Match found
            }
        }
    }

    // If not found in the list, the integral is zero (sparse format)
    return 0.0;
}