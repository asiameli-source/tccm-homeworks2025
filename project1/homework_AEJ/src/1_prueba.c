#include <stdio.h>
#include <trexio.h>
#include <stdlib.h>
#include <stdint.h>

int main() {
  int num = 3;  // Number of atoms
  double coord[][3] = {
    // xyz coordinates in atomic units
    0.  ,  0.        , -0.24962655,
    0.  ,  2.70519714,  1.85136466,
    0.  , -2.70519714,  1.85136466 };
 
struct 
 trexio_exit_code rc;
 
 // Open TREXIO file
 trexio_t* trexio_file = trexio_open("h2o.h5", 'r', TREXIO_AUTO, &rc);
 if (rc != TREXIO_SUCCESS) {
   printf("TREXIO Error: %s\n", trexio_string_of_error(rc));
   exit(1);
 }

// Write coord
 trexio_write_nucleus_coord(trexio_file, &coord[0][0]);
// Reading Nuclear repulsion energy
 double energy;
 rc = trexio_read_nucleus_repulsion(trexio_file, &energy);
 //---check the return code---
 if (rc != TREXIO_SUCCESS) {
  printf("TREXIO Error reading nuclear repulsion energy: \n%s\n", trexio_string_of_error(rc));
  exit(1);
 }

//
// Obtain Number of up electrons (Nocc)
 int32_t n_up = 0;
 rc = trexio_read_electron_up_num(trexio_file, &n_up);
 if (rc != TREXIO_SUCCESS) {
   printf("TREXIO Error: %s\n", trexio_string_of_error(rc));
   exit(1);
 }
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
///
/// One-electron integras (core Hamiltonian)
 double* hcore =malloc((int_64_t)mo_num * mo_num * sizeof(double));
 rc = trexio_read_mo_1e_int_core_hamiltonian(trexio_file, data);
 if (rc != TREXIO_SUCCESS) {
   printf("TREXIO Error: %s\n", trexio_string_of_error(rc));
   exit(1);
 }
// Reading two-electron integrals
/// Number of non-zero integrals
 rc = trexio_read_mo_2e_int_eri_size(trexio_t* const trexio_file, int64_t* const n_integrals);
 if (rc != TREXIO_SUCCESS) {
   printf("TREXIO Error: %s\n", trexio_string_of_error(rc));
   exit(1);
 }
///
///
/// Allocate memory
 int32_t* const index = malloc(4 * n_integrals * sizeof(int32_t));
 if (index == NULL) {
  fprintf(stderr, "Malloc failed for index");
  exit(1);
 }

 double* const value = malloc(n_integrals * sizeof(double));
 if (value == NULL) {
  fprintF(stderr, "Malloc failed for value");
  exit(1);
 }
///
///
/// Read the integrals from the file
 rc = trexio_read_mo_2e_int_eri(trexio_t* const file, const int64_t offset_file, int64_t* const buffer_size, int32_t* const index, double* const value);
 if (rc != TREXIO_SUCCESS) {
   printf("TREXIO Error: %s\n", trexio_string_of_error(rc));
   exit(1);
 }
///
 offset_file=0
 buffer_size = n_integrals

 int i = index [4*n+0];
 int j = index [4*n+1];
 int k = index [4*n+2];
 int l = index [4*n+3];
 double integral = value[n];

///
// Close TREXIO file
 rc = trexio_close(trexio_file);
 if (rc != TREXIO_SUCCESS) {
  printf("TREXIO Error: %s\n", trexio_string_of_error(rc));
  exit(1);
 }
 trexio_file = NULL;
}
