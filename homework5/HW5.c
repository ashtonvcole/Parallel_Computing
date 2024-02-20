/*
 * HW4.c
 * Even/odd transposition sort
 * Ashton Cole
 * COE 379L: Parallel Computing
 */

#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>

#include "tools.h"

int main(int argc, char **argv) {
	MPI_Init(&argc, &argv);

	// MPI Variables
	MPI_Comm comm = MPI_COMM_WORLD;
	int nprocs, procno;
	MPI_Comm_size(comm, &nprocs);
	MPI_Comm_rank(comm, &procno);

	int N = 100;
	double indata[N], outdata[N];
	for (int i = 0; i < N; i++) indata[i] = 1;

	double leftdata = 0., rightdata = 0.;
	int send_to, recv_fr;
	MPI_Request requests[4];

	// Get from left, send rightmost to right
	send_to = procno < nprocs - 1 ? procno + 1 : MPI_PROC_NULL;
	recv_fr = procno > 0 ? procno - 1 : MPI_PROC_NULL;

	MPI_Isend(
			&(indata[N - 1]),
			1,
			MPI_DOUBLE,
			send_to,
			0,
			comm,
			&(requests[0])
		 );

	MPI_Irecv(
			&leftdata,
			1,
			MPI_DOUBLE,
			recv_fr,
			0,
			comm,
			&(requests[1])
		 );
	
	// Get from right, send leftmost to left
	send_to = procno > 0 ? procno - 1 : MPI_PROC_NULL;
	recv_fr = procno < nprocs - 1 ? procno + 1 : MPI_PROC_NULL;

	MPI_Isend(
			&(indata[0]),
			1,
			MPI_DOUBLE,
			send_to,
			0,
			comm,
			&(requests[2])
		 );

	MPI_Irecv(
			&rightdata,
			1,
			MPI_DOUBLE,
			recv_fr,
			0,
			comm,
			&(requests[3])
		 );

	// Wait for all communications to finish
	MPI_Wait(&(requests[0]), MPI_STATUS_IGNORE);
	MPI_Wait(&(requests[1]), MPI_STATUS_IGNORE);
	MPI_Wait(&(requests[2]), MPI_STATUS_IGNORE);
	MPI_Wait(&(requests[3]), MPI_STATUS_IGNORE);

	// Summation
	for (int i = 0; i < N; i++) {
		if (i == 0) {
			outdata[i] = leftdata + indata[i] + indata[i + 1];
		} else if (i == N - 1) {
			outdata[i] = indata[i - 1] + indata[i] + rightdata;
		} else {
			outdata[i] = indata[i - 1] + indata[i] + indata[i + 1];
		}
	}

	// Check results
	double answer[N];
	for (int i = 0; i < N; i++) {
		if ((procno == 0 && i == 0) ||
		    (procno == nprocs - 1 && i == N - 1)) {
			answer[i] = 2.;
		} else {
			answer[i] = 3.;
		}
	}

	int error_test = array_error(answer, outdata, N);
	print_final_result(error_test, comm);

	MPI_Finalize();

	return 0;
}
