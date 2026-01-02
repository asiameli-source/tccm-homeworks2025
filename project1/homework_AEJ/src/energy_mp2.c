#include <stdio.h>
#include <stdlib.h> //due to the exit(1)
#include <stdint.h> //due to the int32_t
#include <trexio.h>
#include<math.h>

// function to swap indices
static inline void swap_int(int *a, int *b) {
  int t = *a;
  *a = *b;
  *b = t;
}

// index mapping helper
// turn a 4-index object into one array address
static inline size_t idx_ijab(int i,int j,int a,int b,int n_occ,int n_virt){
  return (size_t)b + (size_t)n_virt * ((size_t)a + (size_t)n_virt * ((size_t)j + (size_t)n_occ*(size_t)i));
}

// two-electron integrals obey 8-fold permutational symmetry
// function to handle the symmetry
// put the integrals all in the same order
static void canonicalize_eri(int *p, int *q, int *r, int *s) { // 'void' to rewrite indices
        // sort each pair
        if (*q < *p) swap_int(p,q);
        if (*s < *r) swap_int(r,s);
        // compare the two pairs (occupied and virtual)
        if ((*r < *p) || (*r == *p && *s < *q)) {
                swap_int(p,r);
		swap_int(q,s);
        }
}

int main(){
	// open a TREXIO file
	trexio_exit_code rc;
	trexio_t* trexio_file = trexio_open("h2o.h5", 'r', TREXIO_AUTO, &rc);
	if(rc != TREXIO_SUCCESS){
		printf("TREXIO Error: %s\n", trexio_string_of_error(rc));
		exit(1);
	}

	int32_t n_up = -1; // physically impossible for an electron count
			   // know immediately the value was or not set correctly

	// obtaining the number of occupied orbitals for closed-shell system
	//trexio_exit_code trexio_read_electron_up_num(trexio_t* const trexio_file, int32_t* const n_up);
	rc = trexio_read_electron_up_num(trexio_file, &n_up);
	if (rc != TREXIO_SUCCESS) {
                printf("TREXIO Error(read electron.up_num): %s\n", trexio_string_of_error(rc));
                exit(1);
        }
	printf("n_up = %d\n", n_up);
	
	int32_t n_occ = n_up;
	printf("n_occ for a closed shell system= %d\n", n_occ);

	//read the number of molecular orbitals
	// trexio_exit_code trexio_read_mo_num(trexio_t* const trexio_file, int32_t* const mo_num);
	int32_t mo_num = -1;
	rc = trexio_read_mo_num(trexio_file, &mo_num);
	if (rc != TREXIO_SUCCESS) {
                printf("TREXIO Error(read mo.num): %s\n", trexio_string_of_error(rc));
                exit(1);
        }
	printf("mo_num = %d\n", mo_num);

	//obtaining the number of virtual orbitals for closed-shell system
	int32_t n_virt = mo_num - n_occ;
	printf("n_virt = %d\n", n_virt);

	//read orbital energies
	//trexio_exit_code trexio_read_mo_energy(trexio_t* const file, double* const mo_energy);
	double* mo_energy = (double*) malloc((size_t)mo_num * sizeof(double));
	rc = trexio_read_mo_energy(trexio_file, mo_energy);
	if (rc != TREXIO_SUCCESS) {
                printf("TREXIO Error(read mo.energy): %s\n", trexio_string_of_error(rc));
                exit(1);
        }
	printf("mo_energies: \n");
	for (int i = 0; i < mo_num; i++) {
		printf(" eps[%d] = %.10f\n", i, mo_energy[i]);
	}

	// read two-electron integrals
	// first, get the number of non-zero integrals
	// trexio_exit_code trexio_read_mo_2e_int_eri_size(trexio_t* const trexio_file, int64_t* const n_integrals);
	int64_t n_integrals = -1;
	rc = trexio_read_mo_2e_int_eri_size(trexio_file, &n_integrals);
	if (rc != TREXIO_SUCCESS) {
		printf("TREXIO Error(read n_integrals): %s\n", trexio_string_of_error(rc)); 
	}
	printf("n_integrals = %ld\n", n_integrals);
	// allocate memory for indices and values
	int32_t* index = malloc(4 * n_integrals * sizeof(int32_t));
	if (index == NULL) {
		fprintf(stderr, "Malloc failed for index");
		exit(1);
	}
	double* value = malloc(n_integrals * sizeof(double));
	if (value == NULL) {
		fprintf(stderr, "Malloc failed for value");
		exit(1);
	}
	// then, read two electron integrals from the file
	int64_t offset_file = 0;
	int64_t buffer_size = n_integrals;
	rc = trexio_read_mo_2e_int_eri(trexio_file, offset_file, &buffer_size, index, value);
	if (rc != TREXIO_SUCCESS) {
		printf("TREXIO Error(read eri): %s\n", trexio_string_of_error(rc));
	}

	// store 4D tensor as 1D array, object G(i, j, a, b) of the four indexes
	size_t G_size = (size_t)n_occ * n_occ * n_virt * n_virt; // compute double numbers to store
	double* G = calloc(G_size, sizeof(double)); // allocate memory for the MP2 integral
        if (G == NULL) { fprintf(stderr, "calloc failed for G\n"); exit(1); }

	for (int64_t n = 0; n < buffer_size; n++) {
		int i = index[4*n + 0];
		int j = index[4*n + 1];
		int k = index[4*n + 2];
		int l = index[4*n + 3];
		double eri = value[n];
		//if (n < 40) {
	       	//	printf("ERI %ld: (%d %d | %d %d) = %.12e\n", (long)n, i, j, k, l, eri);
		//}
		canonicalize_eri(&i, &j, &k, &l);
		static int printed = 0;
		// keep only occ-occ | virt-virt
		if (i < n_occ && j < n_occ && k >= n_occ && l >= n_occ) {
			//printf("<%d %d | %d %d> = %.12e\n", i,j,k,l,eri);
			int a = k - n_occ;
			int b = l - n_occ;
			size_t pos = idx_ijab(i, j, a, b, n_occ, n_virt); // address calculation
			G[pos] += eri; // if the same integral arrives twice, you accumulate it rather than overwrite
		        double ijab = G[idx_ijab(i,j,a,b,n_occ,n_virt)];
		}	
	}

	// close the TREXIO file
	rc = trexio_close(trexio_file);
	if (rc != TREXIO_SUCCESS) {
		printf("TREXIO Error: %s\n", trexio_string_of_error(rc));
		exit(1);
	}
	trexio_file = NULL;
	// where each TREXIO function returns an exit code, indicating success

	// compute the energy_MP2 equation
	
	// convenient access macro: a,b are 0..n_virt-1 (virtual indices)
	#define G4(i,j,a,b)  ( G[idx_ijab((i),(j),(a),(b), n_occ, n_virt)] )
	double energy_mp2 = 0;
	for (int i=0; i < n_occ; i++){
		for (int j=0; j < n_occ; j++){
			for (int a=0; a < n_virt; a++){
				for (int b=0; b < n_virt; b++){
					double ijab = G4(i, j, a, b);
					double ijba = G4(i, j, b, a);
					double denom = mo_energy[i] + mo_energy[j] - mo_energy[a + n_occ] - mo_energy[b + n_occ];
					energy_mp2 += ijab * (2.0*ijab - ijba) / denom;
				}
			}
		}
	}
	printf("Energy_MP2 = %.12f\n", energy_mp2);
	// prevent memory leaks
	free(mo_energy);
        free(G);
        free(index);
        free(value);
	return 0;
}

