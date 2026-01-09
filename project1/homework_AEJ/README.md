# Project #1 — Hartree–Fock and MP2 Energies in C using TREXIO

This project implements the computation of the closed-shell Hartree–Fock (HF) total energy and the MP2 correlation energy in the molecular-orbital (MO) basis. Input data (MO one- and two-electron integrals, MO energies, electron counts, nuclear repulsion energy) are read from TREXIO files (`.h5`).

The code is written in C and links against the TREXIO C library.

Project context: Advanced Computational Techniques / TCCM homework, Project #1.

---

## Directory structure

- `src/`
  - `main.c`  
  Main program. Reads data from a TREXIO file, computes the Hartree–Fock
  energy and the MP2 correlation energy, and prints the results.
  - `MP2_energy.c`  
  External function that handles the symmetry of the two-electron integrals during the computations.
  - `MP2_energy.h`
  Header file defining the helper interface for ERI retrieval from TREXIO sparse arrays (index + value).
  Declares utility routines to handle ERI index permutations (8-fold symmetry) and to query specific integrals \f$\langle pq|rs\rangle\f$ needed in HF/MP2. 
  - `makefile`  
  Compilation rules and linking against the TREXIO library.
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

### Hartree–Fock total energy
The Hartree–Fock total energy is computed as:

$$
E_{\mathrm{HF}} = E_{\mathrm{NN}}
+ 2 \sum_{i \in \mathrm{occ}} \langle i|h|i\rangle
+ \sum_{i \in \mathrm{occ}} \sum_{j \in \mathrm{occ}}
\left[ 2\langle ij|ij\rangle - \langle ij|ji\rangle \right]
$$

where:
- \f$E_{\mathrm{NN}}\f$ is the nuclear repulsion energy,
- \f$\langle i|h|i\rangle\f$ are matrix elements of the one-electron Hamiltonian, with \f$h = T + V_{\mathrm{ne}}\f$,
- \f$\langle ij|kl\rangle\f$ are two-electron Coulomb integrals in the MO basis,
- indices \f$i,j\f$ run over occupied spatial orbitals (\f$i,j \in \mathrm{occ}\f$).

Implementation note: in the code, occupied indices `i,j` correspond to MO indices `[0, n_occ-1]` (0-based indexing).

### MP2 energy expression
The MP2 correlation energy is computed as

$$
E_{\mathrm{MP2}} =
\sum_{i,j \in \mathrm{occ}}
\sum_{a,b \in \mathrm{virt}}
\frac{\langle ij|ab\rangle \left[2\langle ij|ab\rangle - \langle ij|ba\rangle\right]}
{\varepsilon_i + \varepsilon_j - \varepsilon_a - \varepsilon_b}
$$

where:
- \f$\varepsilon_p\f$ are MO orbital energies,
- \f$\langle ij|ab\rangle\f$ are two-electron repulsion integrals (ERIs) in the MO basis,
- indices \f$i,j\f$ run over occupied orbitals (\f$i,j \in \mathrm{occ}\f$) and \f$a,b\f$ run over virtual orbitals (\f$a,b \in \mathrm{virt}\f$).

Important implementation note: in the code, occupied indices `i,j` run over `[0, n_occ-1]`, while virtual indices `a,b` are stored as **shifted** indices `[0, n_virt-1]` corresponding to MO indices `[n_occ, mo_num-1]`.

---

## Implementation overview

1.  Open TREXIO file in read mode
2.  Read the nuclear repulsion energy \(E_{NN}\)
3.  Determine the occupied space
    - read the number of spin-up electrons;
4.  Read the MO dimension and define the virtual space
    - MO numbers
    - compute number of virtual orbitals
5.  Read MO orbital energies
    - allocate
6.  Read one-electron integrals (core Hamiltonian, MO basis)
    - allocate the dense matrix
7.  Read sparse two-electron integrals (ERIs)
    - read the number of stored ERIs
    - allocate sparse buffers
    - read the ERI block
8.  ERI access via 'get_eri(...)' function
9.  Compute the closed-shell Hartree-Fock total energy
10. Compute the MP2 correlation energy
11. Print, close and free
    - print 'E_hf', 'E_MP2' and 'E_hf + E_MP2'
    - Close TREXIO file
    - Free heap allocations

## Quick usage

Typical execution is:

- compile the code using `make` (see `INSTALL.md`)
- run by providing a dataset path, e.g.
  - `../data/h2o.h5`
  - `../data/ch4.h5`
  - etc ...

If your program expects a specific input file name (default is usually `h2o.h5`), ensure this file is present in the right directory where you execute the binary.

```bash
# Example if the file is hardcoded
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
given by README.org with the expected energies.

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
E_HF    = -76.0267987082
Energy_MP2 = -0.203959974098
E_total (HF+MP2) = -76.230758682348
```
## Authors
Asia Meli, Yiyi Yang and Eloá Abreu. 

Project for Advanced Computational Techniques, 2025.

## References 
Szabo, A., & Ostlund, N. S. (1982). Modern quantum chemistry : introduction to advanced electronic structure theory. 

Trex-CoE. (n.d.). TREXIO source code documentation. https://trex-coe.github.io/trexio/

Doxygen homepage. (n.d.). https://www.doxygen.nl/

ChaGPT for any doubts, understand C syntax, diagnose compilation/linking/runtime issues, as a support tool to clarify both theoretical and practical espects of the implementation, assist in drafting and structuring the repository documentation 
