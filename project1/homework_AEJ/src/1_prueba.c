#include <stdio.h>
#include <trexio.h>
#include <stdlib.h>
#include <stdint.h>

int main() {
 
 trexio_exit_code rc;
 
 // Open TREXIO file
 trexio_t* trexio_file = trexio_open("h2o.h5", 'r', TREXIO_AUTO, &rc);
 if (rc != TREXIO_SUCCESS) {
   printf("TREXIO Error: %s\n", trexio_string_of_error(rc));
   exit(1);
 }


// Reading Nuclear repulsion energy
 double energy;
 rc = trexio_read_nucleus_repulsion(trexio_file, &energy);
 //---check the return code---
 if (rc != TREXIO_SUCCESS) {
  printf("TREXIO Error reading nuclear repulsion energy: \n%s\n", trexio_string_of_error(rc));
  exit(1);
 }
 printf("E_nn = %.10f\n", energy);
//
// Obtain Number of up electrons (Nocc)
 int32_t n_up = 0;
 rc = trexio_read_electron_up_num(trexio_file, &n_up);
 if (rc != TREXIO_SUCCESS) {
   printf("TREXIO Error: %s\n", trexio_string_of_error(rc));
   exit(1);
 }
 printf("n_up = %d\n", n_up);
//
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
///
/// One-electron integrals (core Hamiltonian)
 double* data = malloc((int64_t)mo_num * mo_num * sizeof(double));
 if (data == NULL) { fprintf(stderr, "Malloc failed for core\n"); exit(1); }

 rc = trexio_read_mo_1e_int_core_hamiltonian(trexio_file, data);
 if (rc != TREXIO_SUCCESS) {
   printf("TREXIO Error: %s\n", trexio_string_of_error(rc));
   exit(1);
 }
// Reading two-electron integrals
/// Number of non-zero integrals
 int64_t n_integrals = 0;
 rc = trexio_read_mo_2e_int_eri_size(trexio_file, &n_integrals);
 if (rc != TREXIO_SUCCESS) {
   printf("TREXIO Error: %s\n", trexio_string_of_error(rc));
   exit(1);
 }
///
///
/// Allocate memory
 int32_t* index = malloc(4 * n_integrals * sizeof(int32_t));
 double* value = malloc(n_integrals * sizeof(double));
 if (index == NULL || value == NULL) {
  fprintf(stderr, "Malloc failed for eri\n");
  exit(1);
 }
///
///
/// Read the integrals from the file
 int64_t offset_file = 0;
 int64_t buffer_size = n_integrals; 
 rc = trexio_read_mo_2e_int_eri(trexio_file, offset_file, &buffer_size,  index, value);
 if (rc != TREXIO_SUCCESS) {
   printf("TREXIO Error: %s\n", trexio_string_of_error(rc));
   exit(1);
 }
 printf("n_integrals read = %ld\n", (long) buffer_size);
///
///
// Close TREXIO file
 rc = trexio_close(trexio_file);
 if (rc != TREXIO_SUCCESS) {
  printf("TREXIO Error: %s\n", trexio_string_of_error(rc));
  exit(1);
 }
 trexio_file = NULL;

 free(data);
 free(index);
 free(value);
return 0;
}
