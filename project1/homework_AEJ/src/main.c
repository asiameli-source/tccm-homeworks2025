#include <stdio.h>
#include <trexio.h>
#include <stdlib.h>
#include <stdint.h>
#include "MP2_energy.h"

//TREXIO PROGRAM
int main() {

 trexio_exit_code rc;

 // Open TREXIO file
 trexio_t* trexio_file = trexio_open("h2o.h5", 'r', TREXIO_AUTO, &rc);
 if (rc != TREXIO_SUCCESS) {
   printf("TREXIO Error: %s\n", trexio_string_of_error(rc));
   exit(1);
 }


// Reading Nuclear repulsion energy
 double Enn;
 rc = trexio_read_nucleus_repulsion(trexio_file, &Enn);
 //---check the return code---
 if (rc != TREXIO_SUCCESS) {
  printf("TREXIO Error reading nuclear repulsion energy: \n%s\n", trexio_string_of_error(rc));
  exit(1);
 }
 printf("E_nn = %.10f\n", Enn);

// Obtain Number of up electrons (Nocc)
 int32_t n_up = 0;
 rc = trexio_read_electron_up_num(trexio_file, &n_up);
 if (rc != TREXIO_SUCCESS) {
   printf("TREXIO Error: %s\n", trexio_string_of_error(rc));
   exit(1);
 }
 printf("n_up = %d\n", n_up);
// For closed-shell: Nocc = N_up
 int32_t n_occ = n_up;

// Reading one-electron integrals
//
/// Molecular orbital number
 int32_t mo_num = 0;
 rc = trexio_read_mo_num(trexio_file, &mo_num);
 if (rc != TREXIO_SUCCESS) {
   printf("TREXIO Error: %s\n", trexio_string_of_error(rc));
   exit(1);
 }
 printf("mo_num = %d\n", mo_num);

//obtaining the number of virtual orbitals for closed-shell system
 int32_t n_virt = mo_num - n_occ;
 printf("n_virt = %d\n", n_virt);

//read orbital energies
//trexio_exit_code trexio_read_mo_energy(trexio_t* const file, double* const mo_energy);
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
 double* data = malloc((int64_t)mo_num * mo_num * sizeof(double));
 if (data == NULL) { fprintf(stderr, "Malloc failed for core\n"); exit(1); }

 rc = trexio_read_mo_1e_int_core_hamiltonian(trexio_file, data);
 if (rc != TREXIO_SUCCESS) {
   printf("TREXIO Error: %s\n", trexio_string_of_error(rc));
   exit(1);
 }
// Reading two-electron integrals
//
/// Number of non-zero integrals
 int64_t n_integrals = 0;
 rc = trexio_read_mo_2e_int_eri_size(trexio_file, &n_integrals);
 if (rc != TREXIO_SUCCESS) {
   printf("TREXIO Error: %s\n", trexio_string_of_error(rc));
   exit(1);
 }
/// Allocate memory
 int32_t* index = malloc(4 * n_integrals * sizeof(int32_t));
 double* value = malloc(n_integrals * sizeof(double));
 if (index == NULL || value == NULL) {
  fprintf(stderr, "Malloc failed for eri\n");
  exit(1);
 }
/// Read the integrals from the file
 int64_t offset_file = 0;
 int64_t buffer_size = n_integrals;
 rc = trexio_read_mo_2e_int_eri(trexio_file, offset_file, &buffer_size,  index, value);
 if (rc != TREXIO_SUCCESS) {
   printf("TREXIO Error: %s\n", trexio_string_of_error(rc));
   exit(1);
 }
 printf("n_integrals read = %ld\n", (long) buffer_size);

// Compute Hartree-Fock energy
double E1 = 0;
double E2 = 0;
/// Term first sum
for (int32_t i = 0; i < n_occ; i++) {
 double hii = data[(int64_t)i * mo_num + i];
 E1 += 2.0 * hii;
}
/// Term second sum
 for (int32_t i=0; i<n_occ; i++) {
  for (int32_t j=0; j<n_occ; j++) {
   double J = get_eri(i,j,i,j, buffer_size, index, value);
   double K = get_eri(i,j,j,i,buffer_size, index, value);
   E2 += (2.0 * J - K);
 }
}
// Total sum of Hartree-Fock energy
double E_hf = Enn + E1 + E2;



// Compute MP2 correction

// convenient access macro: a,b are 0..n_virt-1 (virtual indices)
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
printf("E_HF    = %.10f\n", E_hf);
printf("Energy_MP2 = %.12f\n", energy_mp2);
printf("E_total (HF+MP2) = %.12f\n", E_hf + energy_mp2);

// Checking energies of the molecular orbitals
// printf("Check eps: eps[n_occ-1]=%.8f eps[n_occ]=%.8f\n", mo_energy[n_occ-1], mo_energy[n_occ]);
// Close TREXIO file
 rc = trexio_close(trexio_file);
 if (rc != TREXIO_SUCCESS) {
  printf("TREXIO Error: %s\n", trexio_string_of_error(rc));
  exit(1);
 }
 trexio_file = NULL;
// Free memory
 free(mo_energy);
 free(data);
 free(index);
 free(value);
 return 0;
}

