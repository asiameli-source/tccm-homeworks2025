# Installation and Compilation Instructions

This project computes the Hartree-Fock energy and MP2 correlation energy using the TREXIO library.

## 1. Prerequisites

To compile and run this program, you must install the HDF5 library and the TREXIO library.

### 1.1. HDF5 Library
TREXIO depends on HDF5. Install it via your package manager:

* **Ubuntu/Debian:**
  ```bash
  sudo apt install libhdf5-dev
* **MACOS:**
brew install hdf5

### 1.2. TREXIO Library (v2.5.0)

You need to manually install TREXIO. Follow these steps:
1. Download the source code: https://github.com/TREX-CoE/trexio/releases/download/v2.5.0/trexio-2.5.0.tar.gz
2. Unzip and install :
   ``` bash
   tar -zxvf trexio-2.5.0.tar.gz
   cd trexio-2.5.0
   ./configure
   make
   sudo make install

(This usually installs the library to /usr/local/lib and headers to /usr/local/include.)

## 2. Compilation

Run the following command on the working directory:
```bash
 gcc -I/usr/local/include -L/usr/local/lib src/*.c -o mp2 -ltrexio -lm
 ```
### 2.1 Flags explanation:

-I/usr/local/include: Look for header files (trexio.h) here.

-L/usr/local/lib: Look for library files here.

src/*.c: Compiles all C source files found in the src folder.

-ltrexio: Links the TREXIO library.

-lm: Links the standard math library.

-o mp2: Creates the executable file named mp2.

## 3. Use of program

Execute the generated mp2 binary file in the shell
```bash
./mp2
```

### Correction of errors

If you receive an error saying libtrexio.so cannot be found, add the library path to your environment variables before running:

```bash
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/lib
./mp2
 ```

