# Homework 4

This assignment implements a parallelized "bubble sort" algorithm. Even and odd processors share and compare values with their neighbors, deciding whether to swap or not, in an alternating fashion. While this is not necessarily the best sorting algorithm, it lends itself well to parallelization, and demonstrates non-blocking message passing using `MPI_Sendrecv()`.

## Setup and Execution

This C code can be compiled using any MPI compiler script or linking the appropriate libraries. To build the project, call `cmake` followed by a path to the source directory. To run, call the appropriate MPI run command, e.g. `mpirun -n 4 <executable_name>`.

## Algorithm

Each processor holds a single item in the list. Again, this lends itself well to parallelization. First even processors and their right neighbor share, compare, and swap values, and then odd processors do the same with their right neighbor. Once two rounds have passed without swapping, the algorithm finishes. The following example shows how this works on 5 processors.

| Round | 0 | 1 | 2 | 3 | Notes |
| --- | --- | --- | --- | --- | --- |
| X | 8 | 3 | 1 | 2 | Start |
| 1 | 3 | 8 | 1 | 2 | Even-odd |
| 2 | 3 | 1 | 8 | 2 | Odd-even |
| 3 | 1 | 3 | 2 | 8 | ... |
| 4 | 1 | 2 | 3 | 8 | ... |
| 4 | 1 | 2 | 3 | 8 | First round without swapping |
| 4 | 1 | 2 | 3 | 8 | Second round without swapping |

The waiting two rounds seems redundant, but it is necessary in a worst case like this.

| Round | 0 | 1 | 2 | 3 | Notes |
| --- | --- | --- | --- | --- | --- |
| X | 3 | 8 | 1 | 2 | Start |
| 1 | 3 | 8 | 1 | 2 | First round without swapping |
