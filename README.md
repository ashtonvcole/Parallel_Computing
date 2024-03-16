# Parallel Computing

This repository contains assignments from COE 379L: Parallel Computing with Victor Eijkhout ([Email](mailto:eijkhout@tacc.utexas.edu), [TACC](https://tacc.utexas.edu/about/staff-directory/victor-eijkhout/), [GitHub](https://github.com/VictorEijkhout)).

## What is parallel computing?

*Parallel computing* refers to writing programs which are solved across multiple concurrent processes. This means that large problems (e.g. a fluid simulation) can be broken into smaller pieces (e.g. zones of the fluid domain) solved separately and more manageably.

## Contents

- [HW1](/HW1): A simple parallel code where individual processes output their rank.
- [HW2](/HW2): A random matrix equation solver using the full Gauss-Jordan method.
- [HW3](/HW3): A program demonstrating `MPI_Scan()`- and `MPI_Gather()`-type routines.
- [HW4](/HW4): A program implementing a bubble sort using the non-blocking send routine `MPI_sendrecv()`.
- [HW5](/homework5): A program conducting a [stencil operation](https://en.wikipedia.org/wiki/Iterative_Stencil_Loops) using the non-blocking send routines `MPI_Isend()` and `MPI_Irecv()`.
- [HW6](/homework6): A program which distributes array elements amongst processors in an interleaved fashion.
- [HW7](/homework7): A program which solves a linear system of equations using the Jacobi method.
