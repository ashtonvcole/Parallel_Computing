/*
 * HW1b.c
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
	int b = 0; // Number of processesi
	char fname[15];
	
	MPI_Comm_rank(MPI_COMM_WORLD, &a);
	MPI_Comm_size(MPI_COMM_WORLD, &b);

	sprintf(fname, "processor%d.txt", a);
	FILE *f = fopen(fname, "w");
	fprintf(f, "Hello from process %d of %d!\n", a, b);
	fclose(f);

	MPI_Finalize();

	return 0;
}
