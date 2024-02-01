/*
 * HW2.c
 * Gauss-Jordan elimination in parallel
 * Homework 2
 * Ashton Cole
 * COE 379L: Parallel Computing
 */

#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
#include <math.h>
#include <unistd.h>

#include "tools.h"

// double array_error(double ref[], double val[], int array_size);
// void print_final_result(int cond, MPI_Comm comm);

int main(int argc, char **argv) {
	// Adapted from Dr. Eijkhout's Skeleton Program
	MPI_Init(&argc, &argv);

	// MPI variables
	MPI_Comm comm = MPI_COMM_WORLD;
	int nprocs;
	int procno;
	MPI_Comm_size(comm, &nprocs);
	MPI_Comm_rank(comm, &procno);

	// Matrix variables
	int N = nprocs;
	double matcol[N]; // matcol[r] Note that this is the procno-th column
	double solution[N]; // solution[r] Reduced column
	double rhs[N]; // rhs[r] Right-hand side column of matrix equation
	double scalings[N]; // scalings[r] Scalings

	// Set up random seed
	srand((int)(procno * (double) RAND_MAX / nprocs));

	// Create random matrix column on this process
	for (int r = 0; r < N; r++) {
		matcol[r] = rand() / (double) RAND_MAX;
		if (r == procno) {
			matcol[r] += 0.5; // Increase the diagonal for stability
		}
	}

	// Make a random vector from the sum of all processes' columns
	MPI_Allreduce(matcol, rhs, N, MPI_DOUBLE, MPI_SUM, comm);
	
	// Gauss-Jordan algorithm
	for (int k = 0; k < N; k++) {
		// Loop across all pivots
		// If this processor has the pivot column, compute scaling vector
		// i.e. how much the pivot row is scaled to make other rows'
		// entry in the pivot column subtract to 0
		if (k == procno) {
			double pivot = matcol[k];
			for (int r = 0; r < N; r++) {
				scalings[r] = matcol[r] / pivot;
			}
		}

		// These scalings need to be communicated to everyone
		// by the kth process number specifically
		MPI_Bcast(scalings, N, MPI_DOUBLE, k, comm);
		
		// Update matrix and RHS
		for (int r = 0; r < N; r++) {
			if (r == k) continue;
			matcol[r] = matcol[r] - scalings[r] * matcol[k];
			rhs[r] = rhs[r] - scalings[r] * rhs[k];
		}
	}

	// Check that all but diagonal is 0
	for (int r = 0; r < N; r++) {
		if (r == procno) continue;
		else if (fabs(matcol[r]) > 1.e-14) {
			printf("Value %f at (%d, %d) should be 0", matcol[r], r, procno);
		}
	}

	// Solve in parallel
	// Recall that this process controls a column
	// So you can sum this scaled column with those from other processes
	// Note that only the diagonal is nonzero
	double local_solution = rhs[procno] / matcol[procno];
	MPI_Allgather(&local_solution, 1, MPI_DOUBLE, solution, 1, MPI_DOUBLE, comm);

	// Check solution
	double answer[N];
	for (int i = 0; i < N; i++) answer[i] = 1.; // RHS was from a sum
	int error_test = array_error(answer, solution, N);
	print_final_result(error_test, comm);

	MPI_Finalize();

	return 0;
}
