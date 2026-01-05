# Project #1 — Hartree–Fock and MP2 Energies in C using TREXIO

This project implements the computation of the closed-shell Hartree–Fock (HF) total energy and the MP2 correlation energy in the molecular-orbital (MO) basis. Input data (MO one- and two-electron integrals, MO energies, electron counts, nuclear repulsion energy) are read from TREXIO files (`.h5`).

All required quantities (one- and two-electron integrals, orbital energies, nuclear repulsion energy) are read from TREXIO (`.h5`) files.The code is written in C and links against the TREXIO C library.

Project context: Advanced Computational Techniques / TCCM homework, Project #1.

---

## Directory structure

- `src/`
  - `main.c`  
  Main program. Reads data from a TREXIO file, computes the Hartree–Fock
  energy and the MP2 correlation energy, and prints the results.
  - `energy_hf.c` — HF energy evaluation (closed-shell)
  - `energy_mp2.c` — MP2 correlation energy evaluation (closed-shell).
   `makefile`  
  Compilation rules and linking against the TREXIO library.
Note: additional source files present in the repository are not used in the
final build and correspond to previous development or testing stages.
- `data/`
  - `*.h5` — TREXIO datasets provided for testing (e.g., `h2o.h5`, `ch4.h5`, ...)
- `INSTALL.md`
  - compilation and run instructions
- `AUTHORS`
  - contributors list
- `LICENSE`
  - license of this code (e.g., MIT)

---

## Requirements

- C compiler (tested with `gcc`)
- TREXIO installed (headers and shared library)
  - Documentation: https://trex-coe.github.io/trexio/
- Standard C libraries

---

## Theory and conventions

### Closed-shell orbital counting
For a closed-shell reference, the number of occupied **spatial** MOs is equal to the number of spin-up electrons:

- `n_occ = n_up`
- Virtual orbitals: `n_virt = mo_num - n_occ`

### Two-electron integrals and symmetry
TREXIO stores MO two-electron integrals (ERIs) in a sparse representation and exploits the 8-fold permutational symmetry. Therefore, only one representative of each symmetry class is stored in the file. During computation, the code must canonicalize/reorder indices consistently so that equivalent integrals are mapped to a unique location in memory.

### MP2 energy expression (closed-shell, MO basis)
The MP2 correlation energy is computed as

$$
E_{\mathrm{MP2}} =
\sum_{i,j \in \mathrm{occ}}
\sum_{a,b \in \mathrm{virt}}
\frac{(ij|ab)\left[2(ij|ab)-(ij|ba)\right]}
{\varepsilon_i+\varepsilon_j-\varepsilon_a-\varepsilon_b}
$$

where \(\varepsilon_p\) are MO orbital energies and \((ij|ab)\) are ERIs in the MO basis.

Important implementation note: in the code, occupied indices `i,j` run over `[0, n_occ-1]`, while virtual indices `a,b` are stored as **shifted** indices `[0, n_virt-1]` corresponding to MO indices `[n_occ, mo_num-1]`.

---

## Implementation overview

1. Open TREXIO file in read mode.
2. Read scalars:
   - `electron_up_num` → `n_occ`
   - `mo_num` → number of MOs
   - (for HF) `nucleus_repulsion` → \(E_{NN}\)
3. Allocate arrays:
   - `mo_energy[mo_num]`
   - `h_core[mo_num * mo_num]` (if needed for HF)
   - MP2 ERI block `G[n_occ*n_occ*n_virt*n_virt]` (dense 1D array)
4. Read data from TREXIO:
   - MO energies
   - one-electron integrals (HF)
   - sparse ERIs (indices + values)
5. Canonicalize ERI indices and store only the MP2-relevant block `(occ,occ|virt,virt)` into `G`.
6. Close TREXIO (after everything is copied into RAM).
7. Compute HF and/or MP2 energies.
8. Print results and free allocated memory.

## Quick usage

Typical execution is:

- compile the code using `make` or `gcc` (see `INSTALL.md`)
- run by providing a dataset path, e.g.
  - `../data/h2o.h5`
  - `../data/ch4.h5`

If your program expects a specific input file name (default is usually `h2o.h5`), ensure this file is present in the current directory where you execute the binary.

```bash
# Example if the file is hardcoded
cp ../data/h2o.h5 .
make
./MP2_energy
```
## Limitations of the program

- The current implementation assumes a closed-shell reference.
- All MP2-relevant ERIs are stored in memory as a dense array, which limits
  the size of systems that can be treated.
- No parallelization is implemented.

## Validation of the results
The computed HF and MP2 energies were validated against reference values
obtained from standard quantum chemistry packages using the same basis set
and molecular geometry.

## Example output
After following the instructions in `install.md` and running the programme, the following will appear on your screen:
```bash
E_nn = 9.1949655588
n_up = 5
mo_num = 24
n_virt = 19
mo_energies: 
 eps[0] = -20.5504129056
 eps[1] = -1.3367079518
 eps[2] = -0.6993360012
 eps[3] = -0.5665672596
 eps[4] = -0.4931469674
 eps[5] = 0.1855792510
 eps[6] = 0.2562590031
 eps[7] = 0.7893769262
 eps[8] = 0.8543469000
 eps[9] = 1.1634990357
 eps[10] = 1.2003880311
 eps[11] = 1.2532917356
 eps[12] = 1.4446530136
 eps[13] = 1.4762517905
 eps[14] = 1.6747291515
 eps[15] = 1.8673061248
 eps[16] = 1.9349291935
 eps[17] = 2.4530526921
 eps[18] = 2.4905194013
 eps[19] = 3.2856783451
 eps[20] = 3.3390040631
 eps[21] = 3.5105917528
 eps[22] = 3.8660278118
 eps[23] = 4.1475335982
n_integrals read = 13458
E_HF    = -90.5031262233
Energy_MP2 = -0.162000251659
E_total (HF+MP2) = -90.665126474953
```
## Authors
Asia Meli
Yiyi Yang
Eloá Abreu

Project for Advanced Computational Techniques, 2025.

## References 
Szabo, A., & Ostlund, N. S. (1982). Modern quantum chemistry : introduction to advanced electronic structure theory. 

Trex-CoE. (n.d.). TREXIO source code documentation. https://trex-coe.github.io/trexio/

The HDF Group - ensuring long-term access and usability of HDF data and supporting users of HDF technologies. (2024, October 16). The HDF5® Library & File Format - The HDF Group - ensuring long-term access and usability of HDF data and supporting users of HDF technologies. The HDF Group - Ensuring Long-term Access and Usability of HDF Data and Supporting Users of HDF Technologies. https://www.hdfgroup.org/solutions/hdf5/

Doxygen homepage. (n.d.). https://www.doxygen.nl/


