# Installation and Compilation Instructions

This project computes the Hartree-Fock energy and MP2 correlation energy using the TREXIO library.

## 1. Prerequisites

To compile and run this program, you must install the HDF5 library and the TREXIO library.

### 1.1. HDF5 Library
TREXIO depends on HDF5. Install it via your package manager:

* **Ubuntu/Debian:**
  ```bash
  sudo apt install libhdf5-dev
  ```
* **MACOS:**
  ```bash
  brew install hdf5
  ```

### 1.2. TREXIO Library (v2.5.0)

You need to manually install TREXIO. Follow these steps:
1. Download the source code: https://github.com/TREX-CoE/trexio/releases/download/v2.5.0/trexio-2.5.0.tar.gz
2. Unzip and install :
   ```bash
   tar -zxvf trexio-2.5.0.tar.gz
   cd trexio-2.5.0
   ./configure
   make
   sudo make install
   ```
   
(This usually installs the library to /usr/local/lib and headers to /usr/local/include.)

## 2. Compilation
Since the project includes a makefile, you do not need to type long compiler commands.

1. Open your terminal in the project directory.

2. Run the make command:
  ```bash
  make
  ```
This will automatically compile main.c and MP2_energy.c, link them against the TREXIO library, and generate the executable named MP2_energy.

If you need to recompile from scratch, you can run make clean followed by make.

This will automatically compile main.c and MP2_energy.c, link them against the TREXIO library, and generate the executable named . 

## 3. Use of program

### 3.1 Input file

The program need a TREXIO input file named h2o.h5 (water molecule) in the same directory where you run the executable.

### 3.1 Running the program

Execute the generated mp2 binary file in the shell

  ```bash
  ./MP2_energy
  ```

## 4. Correction of errors

If you receive an error saying libtrexio.so cannot be found, add the library path to your environment variables before running:

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


