#include <stdio.h>
#include <trexio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
// function to swap indices
static void swap_int(int *a, int *b) {
  int t = *a;
  *a = *b;
  *b = t;
}
// index mapping helper
// turn a 4-index object into one array
static inline size_t idx_ijab(size_t i,size_t j,size_t a,size_t b,size_t n_occ,size_t n_virt){
  return b + n_virt * a + n_virt * (j + n_occ*i);
}

// two-electron integrals obey 8-fold permutational symmetry
// put the integrals all in the same order
static void canonicalize_eri(int *p, int *q, int *r, int *s) { // 'void' to rewrite indices
/// sort each pair
  if (*q < *p) swap_int(p,q);
  if (*s < *r) swap_int(r,s);
/// compare the two pairs (occupied and virtual)
  if ((*r < *p) || (*r == *p && *s < *q)) {
    swap_int(p,r);
    swap_int(q,s);
  }
}


// Eight permutation symmetry
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
 return 0.0;
}
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

// Object G(i, j, a, b) of the four indexes
  int64_t G_size = n_occ * n_occ * n_virt * n_virt;
  double* G = calloc(G_size, sizeof(double)); // allocate memory for the MP2 integral
  if (G == NULL) { fprintf(stderr, "calloc failed for G\n"); exit(1); }

  for (int64_t n = 0; n < buffer_size; n++) {
    int i = index[4*n + 0];
    int j = index[4*n + 1];
    int k = index[4*n + 2];
    int l = index[4*n + 3];
    double eri = value[n];

    canonicalize_eri(&i, &j, &k, &l);
// keep only occ-occ | virt-virt
  if (i < n_occ && j < n_occ && k >= n_occ && l >= n_occ) {
    int a = k - n_occ;
    int b = l - n_occ;
    size_t pos_ab = idx_ijab(i, j, a, b, n_occ, n_virt); // address calculation
    G[pos_ab] += eri;
  }
    // if the same integral arrives twice, you accumulate it rather than overwrite
//    double ijab = G[idx_ijab(i,j,a,b,n_occ,n_virt)];
 }
 

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
printf("E_HF    = %.10f\n", E_hf);


// Compute MP2 correction

// convenient access macro: a,b are 0..n_virt-1 (virtual indices)
#define G4(i,j,a,b)  ( G[idx_ijab((i),(j),(a),(b), n_occ, n_virt)] )
 double energy_mp2 = 0.0;
  for (int i=0; i < n_occ; i++){
   for (int j=0; j < n_occ; j++){
    for (int a=0; a < n_virt; a++){
     for (int b=0; b < n_virt; b++){
       double ijab = G4(i, j, a, b);
       double ijba = G4(i, j, b, a);
       double denom = mo_energy[i] + mo_energy[j] - mo_energy[a + n_occ] - mo_energy[b + n_occ];

       energy_mp2 +=  (ijab * (2.0*ijab - ijba)) / denom;
       }
      }
     }
    }

printf("Energy_MP2 = %.12f\n", energy_mp2);
printf("E_total (HF+MP2) = %.12f\n", E_hf + energy_mp2);	

// Close TREXIO file
 rc = trexio_close(trexio_file);
 if (rc != TREXIO_SUCCESS) {
  printf("TREXIO Error: %s\n", trexio_string_of_error(rc));
  exit(1);
 }
 trexio_file = NULL;
// Free memory
 free(mo_energy);
 free(G);
 free(data);
 free(index);
 free(value);
 return 0;
}
