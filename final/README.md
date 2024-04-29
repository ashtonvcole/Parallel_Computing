# Linear Acoustic Perturbation Solver

This project implements a two-dimensional finite difference solver, which attempts to model the propagation of acoustic perturbations across a steady, approximately-incompressible flow. The solver is written in C, leveraging the [PETSc](https://petsc.org) linear algebra library. More details about the theory behind the solver can be found in the project report.

## Setup and Execution

### Prerequisites

In addition to a standard C compiler, this project requires an installation of PETSc, and by extension an implementation of MPI. Importantly, the `$PETSC_DIR` and `$PETSC_ARCH` environment variables should be set to the PETSc intallation location and build subdirectory, respectively. 

### Problem Setup

To set up a particular case, constant values are defined in `driver.c`, while functions are defined in `driver_functions.c`. Once these are configured appropriately, the program may be built with CMake.

### Problem Execution

This program may either be executed in serial or parallel. The case output directory will be written on a path relative to the present working directory.

The following command will run the program in serial.

```bash
./acoustic_case
```

The following command will run the program in parallel on 4 processors.

```bash
mpirun -n 4 acoustic_case
```

## Project Structure

## Output Structure


