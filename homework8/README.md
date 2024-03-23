# Homework 8

This assignment demonstrates parallelizing a recursive function algorithm using tasks. It calculates the $n$th item in the Fibonacci sequence. Admittedly, even with parallelization, recursion is not the fastest algorithm.

## Setup and Execution

This C code can be compiled using cmake. To build the project, call `cmake` followed by a path to the source directory. On Mac ARM architecture, the OpenMP libraries may not necessarily be found. With a properly installed GCC compiler, `gcc -o <executable_name> HW8.c -fopenmp` should suffice. To run, set the appropriate OpenMP environment variables and call `./<executable_name> <fibonacci_index>`.
