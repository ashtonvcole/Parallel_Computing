/*
 * HW1a.c
 * Make individual processes print
 * Homework 1
 * Ashton Cole
 * COE 379L: Parallel Computing
 */

#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>

int main(int argc, char **argv) {
	
	MPI_Init(&argc, &argv);

	int a = 0; // Process index
	int b = 0; // Number of processes
	
	MPI_Comm_rank(MPI_COMM_WORLD, &a);
	MPI_Comm_size(MPI_COMM_WORLD, &b);

	printf("Hello from process %d of %d!\n", a, b);

	MPI_Finalize();

	return 0;
}
