/**
 * @file main.c
 * @brief Hartree–Fock (HF) total energy and MP2 correlation energy from a TREXIO file.
 *
 * This program reads all required quantities from a TREXIO dataset (here: "h2o.h5")
 * and computes:
 * - The closed-shell Hartree–Fock (restricted) total energy:
 *   \f[
 *   E_{HF} = E_{nn} + 2\sum_{i}^{occ} h_{ii} + \sum_{ij}^{occ}\left(2(ij|ij) - (ij|ji)\right)
 *   \f]
 * - The canonical MP2 correlation energy:
 *   \f[
 *   E_{MP2} = \sum_{ij}^{occ}\sum_{ab}^{virt} \frac{(ij|ab)\left[2(ij|ab)-(ij|ba)\right]}
 *   {\varepsilon_i + \varepsilon_j - \varepsilon_a - \varepsilon_b}
 *   \f]
 *
 * All one-electron and two-electron integrals are assumed to be expressed in
 * the molecular orbital (MO) basis.
 *
 * The two-electron integrals are stored in the TREXIO file in sparse format
 * (index/value arrays). The retrieval of a specific ERI is delegated to the
 * helper function get_eri(...) declared in "MP2_energy.h".
 *
 * @note The program assumes a closed-shell system where the number of occupied
 *       orbitals is equal to the number of spin-up electrons (n_occ = n_up).
 *
 * @date 2025
 */
#include <stdio.h>
#include <trexio.h>
#include <stdlib.h>
#include <stdint.h>
#include "MP2_energy.h"

//TREXIO PROGRAM
/**
 * @brief Main program entry point.
 *
 * The execution flow is described step-by-step below.
 * @return int Returns 0 if execution is successful.
 */
int main() {

 trexio_exit_code rc;
 /**
  * @brief 1\. Open the TREXIO file.
  * We open the file in read-only mode using the automatic backend detection.
  * @code
  * trexio_t* trexio_file = trexio_open("h2o.h5", 'r', TREXIO_AUTO, &rc);
  * if (rc != TREXIO_SUCCESS) {
  * printf("TREXIO Error: %s\n", trexio_string_of_error(rc));
  * exit(1);
  * }
  * @endcode
  */
 trexio_t* trexio_file = trexio_open("../data/h2o.h5", 'r', TREXIO_AUTO, &rc);
 if (rc != TREXIO_SUCCESS) {
   printf("TREXIO Error: %s\n", trexio_string_of_error(rc));
   exit(1);
 }


/**
  * @brief 2\. Read Nuclear repulsion energy.**
  * This is the classical Coulomb repulsion energy between nuclei ($E_{nn}$).
  * @code
  * double Enn;
  * rc = trexio_read_nucleus_repulsion(trexio_file, &Enn);
  * if (rc != TREXIO_SUCCESS) exit(1);
  * printf("E_nn = %.10f\n", Enn);
  * @endcode
  */
 double Enn;
 rc = trexio_read_nucleus_repulsion(trexio_file, &Enn);
 //---check the return code---
 if (rc != TREXIO_SUCCESS) {
  printf("TREXIO Error reading nuclear repulsion energy: \n%s\n", trexio_string_of_error(rc));
  exit(1);
 }
 printf("E_nn = %.10f\n", Enn);

// Obtain Number of up electrons (Nocc)
/**
  * @brief 3\. Define Electron Occupancy.
  * For a restricted closed-shell system, we read the up-spin electrons and assume $N_{occ} = N_{up}$.
  * @code
  * int32_t n_up = 0;
  * trexio_read_electron_up_num(trexio_file, &n_up);
  * int32_t n_occ = n_up;
  * @endcode
  */
 int32_t n_up = 0;
 rc = trexio_read_electron_up_num(trexio_file, &n_up);
 if (rc != TREXIO_SUCCESS) {
   printf("TREXIO Error: %s\n", trexio_string_of_error(rc));
   exit(1);
 }
 printf("n_up = %d\n", n_up);
// For closed-shell: Nocc = N_up
 int32_t n_occ = n_up;

/**
 * @brief 4\. Read Molecular Orbital (MO) dimensions.
 * We need the total number of MOs to calculate the number of virtual orbitals ($N_{virt} = N_{mo} - N_{occ}$).
 * @code
 * int32_t mo_num = 0;
 * trexio_read_mo_num(trexio_file, &mo_num);
 * int32_t n_virt = mo_num - n_occ;
 * @endcode
 */
 int32_t mo_num = 0;
 rc = trexio_read_mo_num(trexio_file, &mo_num);
 if (rc != TREXIO_SUCCESS) {
   printf("TREXIO Error: %s\n", trexio_string_of_error(rc));
   exit(1);
 }
 printf("mo_num = %d\n", mo_num);

 /**
  * @brief Calculated virtual orbitals.
  */
 int32_t n_virt = mo_num - n_occ;
 printf("n_virt = %d\n", n_virt);

//read orbital energies
//trexio_exit_code trexio_read_mo_energy(trexio_t* const file, double* const mo_energy);
/**
  * @brief 5\. Read Orbital Energies ($\epsilon$).
  * Allocates memory and reads the diagonal elements of the Fock matrix in MO basis.
  * @code
  * double* mo_energy = (double*) malloc((size_t)mo_num * sizeof(double));
  * trexio_read_mo_energy(trexio_file, mo_energy);
  * @endcode
  */
 double* mo_energy = (double*) malloc((size_t)mo_num * sizeof(double));
 if (mo_energy == NULL) { fprintf(stderr, "Malloc failed for mo_energy\n"); exit(1); }

 rc = trexio_read_mo_energy(trexio_file, mo_energy);
 if (rc != TREXIO_SUCCESS) {
  printf("TREXIO Error(read mo.energy): %s\n", trexio_string_of_error(rc));
  exit(1);
 }
 // To see the molecular orbital energies
 printf("mo_energies: \n");
 for (int i = 0; i < mo_num; i++) {
   printf(" eps[%d] = %.10f\n", i, mo_energy[i]);
   }


/// One-electron integrals (core Hamiltonian)
/**
  * @brief 6\. Read Core Hamiltonian ($h_{pq}$).
  * Stored as a dense matrix. We allocate $N_{mo} \times N_{mo}$.
  * @code
  * double* data = malloc((int64_t)mo_num * mo_num * sizeof(double));
  * trexio_read_mo_1e_int_core_hamiltonian(trexio_file, data);
  * @endcode
  */
 double* data = malloc((int64_t)mo_num * mo_num * sizeof(double));
 if (data == NULL) { fprintf(stderr, "Malloc failed for core\n"); exit(1); }

 rc = trexio_read_mo_1e_int_core_hamiltonian(trexio_file, data);
 if (rc != TREXIO_SUCCESS) {
   printf("TREXIO Error: %s\n", trexio_string_of_error(rc));
   exit(1);
 }

/**
 * @brief 7\. Prepare Two-Electron Integrals (ERIs).
 * First, we read the number of non-zero integrals to allocate the buffers.
 * Then we allocate `index` (quadruplets) and `value` arrays.
 * @code
 * int64_t n_integrals = 0;
 * trexio_read_mo_2e_int_eri_size(trexio_file, &n_integrals);
 *
 * int32_t* index = malloc(4 * n_integrals * sizeof(int32_t));
 * double* value = malloc(n_integrals * sizeof(double));
 * @endcode
 */
 int64_t n_integrals = 0;
 rc = trexio_read_mo_2e_int_eri_size(trexio_file, &n_integrals);
 if (rc != TREXIO_SUCCESS) {
   printf("TREXIO Error: %s\n", trexio_string_of_error(rc));
   exit(1);
 }
/**
  * @brief Allocating memory for ERIs...
  */
 int32_t* index = malloc(4 * n_integrals * sizeof(int32_t));
 double* value = malloc(n_integrals * sizeof(double));
 if (index == NULL || value == NULL) {
  fprintf(stderr, "Malloc failed for eri\n");
  exit(1);
 }
/**
  * @brief 8\. Read the ERI chunks.
  * @code
  * rc = trexio_read_mo_2e_int_eri(trexio_file, offset_file, &buffer_size,  index, value);
  * @endcode
  */
 int64_t offset_file = 0;
 int64_t buffer_size = n_integrals;
 rc = trexio_read_mo_2e_int_eri(trexio_file, offset_file, &buffer_size,  index, value);
 if (rc != TREXIO_SUCCESS) {
   printf("TREXIO Error: %s\n", trexio_string_of_error(rc));
   exit(1);
 }
 printf("n_integrals read = %ld\n", (long) buffer_size);

/**
 * @brief 9\. Compute Hartree-Fock Energy.
 *
 * **Part A: One-electron term** ($E_1 = 2 \sum h_{ii}$)
 * @code
 * for (int32_t i = 0; i < n_occ; i++) {
 * double hii = data[(int64_t)i * mo_num + i];
 * E1 += 2.0 * hii;
 * }
 * @endcode
 */
double E1 = 0;
double E2 = 0;
/**
 * @brief Computing E1...
 */
for (int32_t i = 0; i < n_occ; i++) {
 double hii = data[(int64_t)i * mo_num + i];
 E1 += 2.0 * hii;
}
/**
 * @brief **Part B: Two-electron term**
 * ($E_2 = \sum (2J - K)$).
 * We use `get_eri` to fetch integrals from the sparse arrays.
 * @code
 * for (int i=0; i<n_occ; i++) {
 * for (int j=0; j<n_occ; j++) {
 * double J = get_eri(i,j,i,j...);
 * double K = get_eri(i,j,j,i...);
 * E2 += (2.0 * J - K);
 * }
 * }
 * @endcode
 */
 for (int32_t i=0; i<n_occ; i++) {
  for (int32_t j=0; j<n_occ; j++) {
   double J = get_eri(i,j,i,j, buffer_size, index, value);
   double K = get_eri(i,j,j,i,buffer_size, index, value);
   E2 += (2.0 * J - K);
 }
}
/**
 * @brief Total HF Sum.
 */
double E_hf = Enn + E1 + E2;

/**
 * @brief 10\. Compute MP2 Correlation Energy.
 *
 * Loops over occupied (i,j) and virtual (a,b) orbitals.
 * \f[ E_{MP2} = \sum \frac{num}{denom} \f]
 * @code
 * for (int i=0; i < n_occ; i++){
 * for (int j=0; j < n_occ; j++){
 * for (int a=n_occ; a < mo_num; a++){
 * for (int b=n_occ; b < mo_num; b++){
 * // ... calculation ...
 * }
 * }
 * }
 * }
 * @endcode
 */
 double energy_mp2 = 0.0;
  for (int i=0; i < n_occ; i++){
   for (int j=0; j < n_occ; j++){
    for (int a=n_occ; a < mo_num; a++){
     for (int b=n_occ; b < mo_num; b++){
       double ijab = get_eri(i, j, a, b, buffer_size, index, value);
       double ijba = get_eri(i, j, b, a, buffer_size, index, value);
       double denom = mo_energy[i] + mo_energy[j] - mo_energy[a] - mo_energy[b];

       energy_mp2 +=  (ijab * (2.0*ijab - ijba)) / denom;
       }
      }
     }
    }
/**
 * @brief 11\. Final Cleanup.
 * Print results and free memory.
 */
printf("E_HF    = %.10f\n", E_hf);
printf("Energy_MP2 = %.12f\n", energy_mp2);
printf("E_total (HF+MP2) = %.12f\n", E_hf + energy_mp2);

// Checking energies of the molecular orbitals
// printf("Check eps: eps[n_occ-1]=%.8f eps[n_occ]=%.8f\n", mo_energy[n_occ-1], mo_energy[n_occ]);
/**
 * @brief Closing file and freeing memory.
 */
 
 rc = trexio_close(trexio_file);
 if (rc != TREXIO_SUCCESS) {
  printf("TREXIO Error: %s\n", trexio_string_of_error(rc));
  exit(1);
 }
 trexio_file = NULL;

 free(mo_energy);
 free(data);
 free(index);
 free(value);
 return 0;
}
