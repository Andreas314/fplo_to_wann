# fplo\_to\_wann

**fplo_to_wann** is a C++ program that converts `+hamdata` files from **FPLO** into `_hr.dat` files compatible with **Wannier90**, and generates a corresponding `POSCAR` file.

## Features

- Converts **FPLO +hamdata** to **Wannier90 _hr.dat**
- Automatically generates a **POSCAR** file
- Supports optional detailed output using the `-H` flag
- Written in modern C++ using the **Eigen** library

### Options

- `-H i j k`  
  Accepts three integers to print hopping terms in the given cell

### Example

Run in the same directory as +hamdata
./fplo\_to\_wann -H 1 2 3

## Installation
git clone https://gitlab.com/yourusername/fplo\_to\_wann.git
cd fplo\_to\_wann
make

## Requirements

- A C++23-compatible compiler (e.g., `g++`, `clang++`)
- [Eigen](https://eigen.tuxfamily.org/) (header-only linear algebra library)
