# Homework 6

This assignment demonstrates how custom data types can be defined within the MPI API. Specifically, it deals with the case where array elements are intended to be split amongst processors in an interleaved fashion. For example, with five processors, every fifth is sent to the fifth processor, the one right before it is sent to the fourth, &c. This is called an "interleaved scatter." It uses a "strided vector," i.e. a vector type of regular blocks of elements separated by regular spacing. In the prior case, this would mean blocks of one element separated by four.

## Setup and Execution

This C code can be compiled using any MPI compiler script or linking the appropriate libraries. To build the project, call `cmake` followed by a path to the source directory. To run, call the appropriate MPI run command, e.g. `mpirun -n 4 <executable_name>`.
