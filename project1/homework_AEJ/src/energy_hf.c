#include <stdio.h>
#include <trexio.h>

int main() {
 int num = 3;
 
 trexio_exit_code rc;

 // Open TREXIO file
 trexio_t* trexio_file = trexio_open(filename, 'r', TREXIO_AUTO, &rc);
 if (rc != TREXIO_SUCCESS) {
   printf("TREXIO Error: %s\n", trexio_string_of_error(rc));
   exit(1);
}

// Reading Nuclear repulsion energy
//
//
// Obtain Number of up electrons (Nocc)
//
//
// Reading one-electron integrals 
//
/// Molecular orbital number
///
///
/// One-electron integras (core Hamiltonian)
//
// Reading two-electron integrals
//
/// Number of non-zero integrals
///
///
///
/// Allocate memory
///
///
/// Read the integrals from the file
///
///
// Close TREXIO file
rc = trexio_close(trexio_file);
if (rc != TREXIO_SUCCESS) {
 printf("TREXIO Error: %s\n", trexio_string_of_error(rc));
 exit(1);
}
trexio_file = NULL;

