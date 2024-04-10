# Homework 10

This assignment demonstrates how various compiler optimizations impact parallel and serial performance of a simple program that performs summation and multiplication operations on vectors.

## Setup and Execution

This C code can be compiled using cmake. To build the project, call `cmake` followed by a path to the source directory. On Mac ARM architecture, the OpenMP libraries may not necessarily be found. With a properly installed GCC compiler, `gcc -o <executable_name> HW10.c -fopenmp` should suffice. To run, set the appropriate OpenMP environment variables and call `./<executable_name> <vector_size>`.
