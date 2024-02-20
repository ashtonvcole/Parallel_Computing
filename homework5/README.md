# Homework 5

This assignment demonstrates how non-blocking sends can be conducted with the `MPI_Isend()` and `MPI_Irecv()` routines. They are used to share values for a list split across multiple processors, which undergoes a summation operation using a [stencil](https://en.wikipedia.org/wiki/Iterative_Stencil_Loops). This is a simple, one-dimensional version of what occurs in simulations where the computational domain is solved in parallel.

## Setup and Execution

This C code can be compiled using any MPI compiler script or linking the appropriate libraries. To build the project, call `cmake` followed by a path to the source directory. To run, call the appropriate MPI run command, e.g. `mpirun -n 4 <executable_name>`.