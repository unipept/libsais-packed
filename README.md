# Direct construction of sparse suffix arrays with Libsais

This repository implements a memory-efficient method for direct sparse suffix array (SSA) construction using the Libsais library. Unlike traditional approaches that first build a full suffix array before sampling, this method directly constructs the SSA through a simple text transformation. By encoding groups of k characters into compact integer representations, the memory usage reduces significantly.

Key Benefits:
* Lower memory usage (no need for a full suffix array)
* Faster SSA construction (directly processes transformed text)
* Optimized for bioinformatics

How It Works:
1) Bit-pack characters into compact integers while preserving lexicographic order.
2) Use Libsais to construct the suffix array of the transformed text.
3) Achieve direct SSA construction, avoiding unnecessary memory overhead.

## Installation
Clone the repository and build the executable:
```
git clone git@github.com:unipept/unipept-libsais.git
cd unipept-libsais
mkdir build && cd build
cmake ..
make
```
This will generate the build_ssa executable.

## Usage
Run the program with the following syntax:
```
./build_ssa -s <sparseness> [-cdu] <input_file> <output_file>
```
### Arguments:
* -s <sparseness>: Defines the sparseness factor (an integer).
* -d: Treats the input file as DNA data (default assumes proteomic data).
* -c: Enables compressed output using bit-packing.
* -u: If enabled, the program will compute the SSA unoptimized, by computing the full SA and subsampling afterwards.
* <input_file>: Path to the input file containing DNA/protein sequences.
* <output_file>: Path where the sparse suffix array will be saved.

### Example
```
./build_ssa -s 3 input.txt output.ssa
```
This command builds an SSA with sparseness factor 3, treats the input as proteomic data, and uses the optimized algorithm.

## Libsais
This tool contains a modified fork of the `libsais` library. The libsais library is a tool for fast linear time suffix array based on induced sorting algorithm described in the following papers: 
* Ge Nong, Sen Zhang, Wai Hong Chan *Two Efficient Algorithms for Linear Suffix Array Construction*, 2009
* Juha Karkkainen, Giovanni Manzini, Simon J. Puglisi *Permuted Longest-Common-Prefix Array*, 2009
* Nataliya Timoshevskaya, Wu-chun Feng *SAIS-OPT: On the characterization and optimization of the SA-IS algorithm for suffix array construction*, 2014
* Jing Yi Xie, Ge Nong, Bin Lao, Wentao Xu *Scalable Suffix Sorting on a Multicore Machine*, 2020

Original Copyright (c) 2021-2024 Ilya Grebnov <ilya.grebnov@gmail.com>

>The libsais is inspired by [libdivsufsort](https://github.com/y-256/libdivsufsort), [sais](https://sites.google.com/site/yuta256/sais) libraries by Yuta Mori and [msufsort](https://github.com/michaelmaniscalco/msufsort) by Michael Maniscalco.

### Changes
* Removed functionality for computing the Burrows-Wheeler transform and longest common prefix array.
* Removed OpenMP acceleration.
* Removed construction of a 32-bit suffix array.
* Added functionality for contstructing a 64-bit suffix array for a 32-bit input.

## License
Unipept-libsais is released under the [Apache License Version 2.0](LICENSE "Apache license")
