/*
 * HW4.c
 * Even/odd transposition sort
 * Ashton Cole
 * COE 379L: Parallel Computing
 */

#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>

int main(int argc, char **argv) {
	MPI_Init(&argc, &argv);

	// MPI Variables
	MPI_Comm comm = MPI_COMM_WORLD;
	int nprocs, procno;
	MPI_Comm_size(comm, &nprocs);
	MPI_Comm_rank(comm, &procno);

	int root = 0;
	int isroot = procno == root;
	int iseven = procno % 2 == 0;

	// Each processor gets an item in the list
	double *list = isroot ? 
		(double *) malloc(nprocs * sizeof(double)) :
		NULL;
	srand((int) (procno * (double) RAND_MAX / nprocs));
	double item = rand() / (double) RAND_MAX; // From 0 to 1
	double new;
	
	// Gather and display to the user
	printf("procno = %d\n - isroot = %d\n - iseven = %d\n", procno, isroot, iseven);/*
	MPI_Gather(
			&item,
			1,
			MPI_DOUBLE,
			&list,
			1,
			MPI_DOUBLE,
			root,
			comm
		  );
	printf("procno = %d\n - isroot = %d\n - iseven = %d\n", procno, isroot, iseven);*/
	/* if (isroot) {
		printf("List before sort: {");
		int i;
		for (i = 0; i < nprocs - 1; i++) {
			//printf("%f, ", list[i]);
		}
		//printf("%f}\n", list[nprocs - 1]);
	} */

	int done = 0;
	int hnswp; // Has this processor not exchanged values this round?

	int sendto, recvfr;

	int count = 0;
	printf("Proc %d starting loop\n", procno);	
	while (!done && count < 100) {
		count = count + 1;

		// Even-odd: 2i and 2i + 1 might swap
		sendto = iseven ? procno + 1 : MPI_PROC_NULL;
		if (sendto < 0 || sendto >= nprocs) sendto = MPI_PROC_NULL;
		recvfr = !iseven ? procno - 1 : MPI_PROC_NULL;
		if (recvfr < 0 || recvfr >= nprocs) recvfr = MPI_PROC_NULL;
		MPI_Sendrecv(
				&item,
				1,
				MPI_DOUBLE,
				sendto,
				0,
				&new,
				1,
				MPI_DOUBLE,
				recvfr,
				0,
				comm,
				MPI_STATUS_IGNORE
			    );

		// Decide whether to swap
		if (iseven) { // Left
			if (new < item) {
				item = new;
				hnswp = 0;
			} else {
				hnswp = 1;
			}
		} else { // Right
			if (new > item) {
				item = new;
				hnswp = 0;
			} else {
				hnswp = 1;
			}
		}

		// If nobody has swapped
		// If proc0 hnswp and proc1 hnswp and ...
		// We are done
		MPI_Allreduce(
				&hnswp,
				&done,
				1,
				MPI_INT,
				MPI_LAND,
				comm
			     );

		// Save some effort
		if (done) break;

		// Odd-even: 2i + 1 and 2i + 2 might swap
		sendto = !iseven ? procno + 1 : MPI_PROC_NULL;
		if (sendto < 0 || sendto >= nprocs) sendto = MPI_PROC_NULL;
		recvfr = iseven ? procno - 1 : MPI_PROC_NULL;
		if (recvfr < 0 || recvfr >= nprocs) recvfr = MPI_PROC_NULL;
		MPI_Sendrecv(
				&item,
				1,
				MPI_DOUBLE,
				sendto,
				0,
				&new,
				1,
				MPI_DOUBLE,
				recvfr,
				0,
				comm,
				MPI_STATUS_IGNORE
			    );

		// Decide whether to swap
		if (!iseven) { // Left
			if (new < item) {
				item = new;
				hnswp = 0;
			} else {
				hnswp = 1;
			}
		} else { // Right
			if (new > item) {
				item = new;
				hnswp = 0;
			} else {
				hnswp = 1;
			}
		}

		// If nobody has swapped
		// If proc0 hnswp and proc1 hnswp and ...
		// We are done
		MPI_Allreduce(
				&hnswp,
				&done,
				1,
				MPI_INT,
				MPI_LAND,
				comm
			     );
	}

	// Gather and display to the user
	MPI_Gather(
			&item,
			1,
			MPI_DOUBLE,
			&list,
			1,
			MPI_DOUBLE,
			root,
			comm
		  );
	if (isroot) {
		printf("List after sort: {");
		for (int i = 0; i < nprocs - 1; i++) {
			//printf("%f, ", list[i]);
		}
		//printf("%f}\n", list[nprocs - 1]);
		free(list);
	}
  
  	MPI_Finalize();
  
  	return 0;
}
