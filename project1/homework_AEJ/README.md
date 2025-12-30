# Project #1 — Hartree–Fock and MP2 Energies in C using TREXIO

This project implements the computation of the closed-shell Hartree–Fock (HF) total energy and the MP2 correlation energy in the molecular-orbital (MO) basis. Input data (MO one- and two-electron integrals, MO energies, electron counts, nuclear repulsion energy) are read from TREXIO files (`.h5`).

The code is written in C and links against the TREXIO C library.

Project context: Advanced Computational Techniques / TCCM homework, Project #1.

---

## Directory structure

- `src/`
  - `energy_hf.c` — HF energy evaluation (closed-shell)
  - `energy_mp2.c` — MP2 correlation energy evaluation (closed-shell)
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

\[
E_{\mathrm{MP2}} =
\sum_{i,j \in \mathrm{occ}}
\sum_{a,b \in \mathrm{virt}}
\frac{(ij|ab)\left[2(ij|ab)-(ij|ba)\right]}
{\varepsilon_i+\varepsilon_j-\varepsilon_a-\varepsilon_b}
\]

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

---

## Quick usage

Typical execution is:

- compile (see `INSTALL.md`)
- run by providing a dataset path, e.g.
  - `../data/h2o.h5`
  - `../data/ch4.h5`

If your pr

