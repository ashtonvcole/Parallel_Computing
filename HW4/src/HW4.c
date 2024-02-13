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
	MPI_Gather(
			&item, // Note that this is not a pointer
			1,
			MPI_DOUBLE,
			list, // Note that this is a pointer
			1,
			MPI_DOUBLE,
			root,
			comm
		  );
	if (isroot) {
		printf("List before sort: {");
		int i;
		for (i = 0; i < nprocs - 1; i++) {
			printf("%f, ", list[i]);
		}
		printf("%f}\n", list[nprocs - 1]);
	}

	int hnswp; // Has this processor not exchanged values this round?
	int all_hnswp; // Has no processor exchanged valuse this round
	int hnswp_cnt = 0; // After two rounds (even-odd and odd-even), we know we are done
			   // Normally one round of no swaps suffices
			   // Except if you start like this
			   // [1, 3, 2, 4] -> hnswp_cnt = 0
			   // Even- odd
			   // [1, 3, 2, 4] -> hnswp_cnt = 1

	int partner; // Note swapping is mutual
	
	while (hnswp_cnt < 2) {
		// Even-odd: 2i and 2i + 1 might swap	
		new = item; // Necessary reset for non-swaps
		partner = iseven ? procno + 1 : procno - 1;
		if (partner < 0 || partner >= nprocs) partner = MPI_PROC_NULL;
		MPI_Sendrecv(
				&item,
				1,
				MPI_DOUBLE,
				partner,
				0,
				&new,
				1,
				MPI_DOUBLE,
				partner,
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
				&all_hnswp,
				1,
				MPI_INT,
				MPI_LAND,
				comm
			     );

		// If at least one has swapped, reset
                hnswp_cnt = all_hnswp ? hnswp_cnt + 1 : 0;

		// Gather and display to the user
                MPI_Gather(
                                &item,
                                1,
                                MPI_DOUBLE,
                                list,
                                1,
                                MPI_DOUBLE,
                                root,
                                comm
                          );
                if (isroot) {
                        printf("List after even-odd: {");
                        for (int i = 0; i < nprocs - 1; i++) {
                                printf("%f, ", list[i]);
                        }
                        printf("%f}\n", list[nprocs - 1]);
                }
		
		// Save some effort
		if (hnswp_cnt >= 2) break;

		// Odd-even: 2i + 1 and 2i + 2 might swap
		new = item; // Necessary reset for non-swaps
		partner = !iseven ? procno + 1 : procno - 1;
		if (partner < 0 || partner >= nprocs) partner = MPI_PROC_NULL;
		MPI_Sendrecv(
				&item,
				1,
				MPI_DOUBLE,
				partner,
				0,
				&new,
				1,
				MPI_DOUBLE,
				partner,
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
				&all_hnswp,
				1,
				MPI_INT,
				MPI_LAND,
				comm
			     );

		// If at least one has swapped, reset
		hnswp_cnt = all_hnswp ? hnswp_cnt + 1 : 0;

		// Gather and display to the user
	        MPI_Gather(
	                        &item,
	                        1,
	                        MPI_DOUBLE,
	                        list,
	                        1,
	                        MPI_DOUBLE,
	                        root,
	                        comm
	                  );
		if (isroot) {
                        printf("List after odd-even: {");
                        for (int i = 0; i < nprocs - 1; i++) {
                                printf("%f, ", list[i]);
                        }
                        printf("%f}\n", list[nprocs - 1]);
                }
	}

	// Gather and display to the user
	MPI_Gather(
			&item,
			1,
			MPI_DOUBLE,
			list,
			1,
			MPI_DOUBLE,
			root,
			comm
		  );
	if (isroot) {
		printf("List after sort: {");
		for (int i = 0; i < nprocs - 1; i++) {
			printf("%f, ", list[i]);
		}
		printf("%f}\n", list[nprocs - 1]);
		free(list);
	}
  
  	MPI_Finalize();
  
  	return 0;
}
