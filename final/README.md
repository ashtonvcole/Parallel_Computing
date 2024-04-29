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

- [`src/`](src/): holds the source code for this project
	- [cmake-modules/](src/cmake-modules/): Ancillary scripts for CMake.
		- [FindPETSc.cmake](src/cmake-modules/FindPetsc.cmake): Passes PETSc environment variables to the main CMake script.
	- [`CMakeLists.txt`](src/CMakeLists.txt)
	- [`driver_functions.c`](src/driver_functions.c): Defines all input functions for a case.
	- [`driver_functions.h`](src/driver_functions.h)
	- [`driver.c`](src/driver.c): Defines all input constants for a case, as well as the main function.
	- [`lib/`](lib/): Holds solver functions.
		- [`CMakeLists.txt`](src/lib/CMakeLists.txt)
		- [`acoustic_problem.c`](src/lib/acoustic_problem.c): Defines all structs and functions needed to solve a case.
		- [`acoustic_problem.h`](src/lib/acoustic_problem.h)
- [`doc/`](doc/): Holds the final report for this project.
	- [`report.tex`](doc/report.tex)
	- [`report.pdf`](doc/report.pdf): The final report for this project, which further explains the theory, implementation, and results of this program.

## Output Structure

To keep your file system clean, every case automatically writes its output to its own new directory.

- `metadata.json`: A file containing information about the domain and output settings.
- `n/`: Integer-numbered folders consisting of the results from the *n*th time step.
	- `t`: A simple data file of the time associated with the time step.
	- `z`: A simple data file of the solution vector associated with the time step. More details on how this vector is organized can be found in the [final report](doc/report.pdf).
