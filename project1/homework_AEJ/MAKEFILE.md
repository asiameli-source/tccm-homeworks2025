# Makefile Configuration

This project uses a standard **Makefile** to automate the compilation of the C source code and its linking against the **TREXIO** chemistry library.

### Project Information
* **Authors:** Asia Meli, Yiyi Yang, Eloá Abreu
* **Date:** 2025
* **Purpose:** Generates the `MP2_energy` executable for Hartree-Fock and MP2 calculations.

---

## 1. Compilation Variables

We define the necessary tools and paths to locate the external libraries.

| Variable | Value | Description |
| :--- | :--- | :--- |
| **CC** | `gcc` | The C compiler to use. |
| **CFLAGS** | `-I/usr/local/include` | Defines where to find the TREXIO header files (`.h`). |
| **LDFLAGS**| `-L/usr/local/lib` | Defines where to find the library binary files (`.so` or `.a`). |
| **LDLIBS** | `-ltrexio` | Instructs the linker to specifically link the `libtrexio` library. |

---

## 2. Source Files and Objects

The system automatically detects what needs to be compiled based on this list:

```makefile
SRC = MP2_energy.c main.c
OBJ = $(patsubst %.c, %.o, $(SRC))
```
Main Target: MP2_energy
This is the final goal. It links all object files (.o) with the defined libraries.

```makefile
MP2_energy: $(OBJ)
    $(CC) -o MP2_energy $(OBJ) $(LDFLAGS) $(LDLIBS)
```
Compiles each source file individually without linking (using the -c flag).

```makefile
%.o: %.c
    $(CC) -c $< $(CFLAGS)
```

Removes generated files to ensure a fresh compilation from scratch.


```makefile
clean:
    rm -f $(OBJ) MP2_energy
```
