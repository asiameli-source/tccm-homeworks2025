# Installation and Compilation Instructions

This project computes the Hartree-Fock energy and MP2 correction using molecular integrals stored in the TREXIO library.

---

## 1. Prerequisites

To compile and run this program, you need:
- **HDF5** with TREXIO dependency
- **TREXIO** (tested with v2.5.0)
- A **C compiler** 


### 1.1. Installing HDF5 Library

Install it via your package manager:

* **Ubuntu/Debian**:
  ```bash
  sudo apt install libhdf5-dev
  ```
* **macOS** (Homebrew):
  ```bash
  brew install hdf5
  ```

### 1.2. Installing TREXIO (v2.5.0)

You need to manually install TREXIO. 
Steps:
1. Download the source code: 
```bash
wget https://github.com/TREX-CoE/trexio/releases/download/v2.5.0/trexio-2.5.0.tar.gz
```
2. Unzip and install :
   ```bash
   tar -zxvf trexio-2.5.0.tar.gz
   cd trexio-2.5.0
   ./configure
   make
   sudo make install
   ```
   
By default, this installs the libraries to `/usr/local/lib` and headers to `/usr/local/include`.

## 2. Compilation

A ``makefile`` is provided in this project,then you do not need to type long compiler commands.

1. Open your terminal in the project directory.

2. Run:
  ```bash
  make
  ```
This compiles ``main.c`` and ``MP2_energy.c``, link them against the TREXIO library, and generate the executable **MP2_energy**.

To rebuild from scratch, run:
 ```bash
 make clean
 make
 ```

## 3. Running the program

### 3.1 Input file

The program expects a TREXIO input file named:
```bash
 molecule.h5
```
In this project, the h2o.h5 (water molecule) was used as an example to compute the energy. The file is added in the same directory where the executable was runned.

### 3.2 Execute

Execute the generated mp2 binary file:

  ```bash
  ./MP2_energy
  ```
The output prints the HF energy, MP2 correlation energy, and the total HF+MP2 energy.

## 4. Correction of errors

### 4.1 libtrexio.so / dynamic library not found (Linux)
If the runtime linker cannot find TREXIO, add the library path to your environment variables before running:

* **Ubuntu/Debian:**
  ```bash
  export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/lib
  ./MP2_energy
  ```
* **MACOS:**
  ```bash
  export DYLD_LIBRARY_PATH=$DYLD_LIBRARY_PATH:/usr/local/lib
  ./MP2_energy
  ```
