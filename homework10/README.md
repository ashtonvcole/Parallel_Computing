# Homework 10

This assignment demonstrates how various compiler optimizations impact parallel and serial performance of a simple program that performs summation and multiplication operations on vectors.

## Setup and Execution

This C code can be compiled using cmake. To build the project, call `cmake` followed by a path to the source directory. On Mac ARM architecture, the OpenMP libraries may not necessarily be found. With a properly installed GCC compiler, `gcc -o <executable_name> HW10.c -fopenmp` should suffice. To run, set the appropriate OpenMP environment variables and call `./<executable_name> <vector_size>`.

### Compiler and Optimization Tests

For this assignment, the Intel 19, Intel 24, and GNU 11 C compilers were used. Their respective commands are `icc`, `icx`, and `gcc`. The optimization flags used in testing were `-O0`, `-O1`, `-O2`, and `-O3`, with a higher number representing additional steps attempting to decrease runtime, at the costs of compilation time and potentially unexpected behavior. The `-qopt-report=3` flag may additionally be used on the Intel compilers to list what optimizations exactly have been used. These reports have been included in the repository for reference.
