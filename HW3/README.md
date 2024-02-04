# Homework 3

This assignment shows how scan and gather routines can be used to assemble an array split across multiple processors. This is relevant in parallel computing, since in some problems, large data structures (e.g. solution nodes on a mesh grid) are split across multiple processors for calculations, and hace to be re-assembled post facto.

## Setup and Execution

This C code can be compiled using any MPI compiler. It is built with the CMake system. To buld the project, call `cmake` followed by a path to the source directory. To run, call the appropriate MPI run command, e.g. `mpirun <executable_name> -n 6`.
