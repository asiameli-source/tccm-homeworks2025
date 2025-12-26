#include <stdio.h>
#include <stdlib.h> //due to the exit(1)
#include <stdint.h> //due to the int32_t
#include <trexio.h>
int main(){
	// open a TREXIO file
	trexio_exit_code rc;
	trexio_t* trexio_file = trexio_open("../data/ch4.h5", 'r', TREXIO_AUTO, &rc);
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
	for (int i = 1; i < mo_num + 1; i++) {
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
	// then, read integrals from the file
	int64_t offset_file = 0;
	int64_t buffer_size = n_integrals;
	rc = trexio_read_mo_2e_int_eri(trexio_file, offset_file, &buffer_size, index, value);
	if (rc != TREXIO_SUCCESS) {
		printf("TREXIO Error(read eri): %s\n", trexio_string_of_error(rc));
	}

	for (int64_t n = 1; n < buffer_size + 1; n++) {
		int i = index[4*n + 0];
		int j = index[4*n + 1];
		int k = index[4*n + 2];
		int l = index[4*n + 3];
		double eri = value[n];
		//if (n < 40) {
	       	//	printf("ERI %ld: (%d %d | %d %d) = %.12e\n", (long)n, i, j, k, l, eri);
		//}
	}

	// close the TREXIO file
	rc = trexio_close(trexio_file);
	if (rc != TREXIO_SUCCESS) {
		printf("TREXIO Error: %s\n", trexio_string_of_error(rc));
		exit(1);
	}
	trexio_file = NULL;
	// where each TREXIO function returns an exit code, indicating success
	return 0;	

}
