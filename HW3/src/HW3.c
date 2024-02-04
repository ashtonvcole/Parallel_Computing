/*
 * HW3.c
 * Scan and Gather Operations
 * Homework 3
 * Ashton Cole
 * COE 379L: Parallel Computing
 */

#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
#include <math.h>
#include <unistd.h>

#include "tools.h"



#define max(a, b) ((a) > (b) ? (a) : (b)) // Macro for maximum of two numbers



int main(int argc, char **argv) {
	// Adapted from Dr. Eijkhout's Skeleton Program
	MPI_Init(&argc, &argv);

	// MPI Variables
	MPI_Comm comm = MPI_COMM_WORLD;
	int nprocs, procno;
	MPI_Comm_size(comm, &nprocs);
	MPI_Comm_rank(comm, &procno);

	// Set a unique random seed
	// Seeds from 0 to RAND_MAX (upper-exclusive), evenly spaced
	srand((int)(procno * (double) RAND_MAX / nprocs));

	// Determine index elements
	int my_number_of_elements = rand() % nprocs; // Random number from 0 to nprocs - 1
	int my_first_index = 0; // Not set yet
	MPI_Exscan( // Scan from 0 to procno - 1
			&my_number_of_elements, // Operands
			&my_first_index, // Where to put result
			1, // Size of operand/result
			MPI_INT, // Type
			MPI_SUM, // Operation
			comm // Communicator
		  );
	printf(
			"Proc %3d has %3d elements, range [%4d,%4d)\n",
			procno,
			my_number_of_elements,
			my_first_index,
			my_first_index + my_number_of_elements
	      );

	// Fill in local array
	int *my_elements = (int *) malloc(max(my_number_of_elements, 1) * sizeof(int)); // Array
	for (int i = 0; i < my_number_of_elements; i++) {
		my_elements[i] = my_first_index + i; // Value is global index
	}

	// Get total number of elements
	int total_number_of_elements;
	MPI_Reduce(
			&my_number_of_elements, // Operands
			&total_number_of_elements, // Where to put result
			1, // Size of operand/result
			MPI_INT, // Type
			MPI_SUM, // Operation
			0, // Process number at which to put the result
			comm // Communicator
		  );
	if (procno == 0) {
		printf("Total number of elements: %d\n", total_number_of_elements);
	}

	// Gather local arrays into a global one at process 0
	// Start by getting an array of the local array sizes
	int *size_buffer = NULL;
	if (procno == 0) size_buffer = (int*) malloc(nprocs * sizeof(int));
	MPI_Gather(
			&my_number_of_elements, // Items to gather
			1, // Size of items per process (uniform)
			MPI_INT, // Type
			size_buffer, // Where to put result
			1, // Size to expect incoming per process (uniform)
			MPI_INT, // Type
			0, // Process number at which to put the result
			comm // Communicator
		  );

	// Offsets/displacements
	int *displ_buffer = NULL;
	if (procno == 0) displ_buffer = (int *) malloc(nprocs * sizeof(int));
	MPI_Gather(
			&my_first_index, // Items to gather
                        1, // Size of items per process (uniform)
                        MPI_INT, // Type
                        displ_buffer, // Where to put result
                        1, // Size to expect incoming per process (uniform)
                        MPI_INT, // Type
                        0, // Process number at which to put the result
                        comm // Communicator
		  );
	
	// Now gather the local arrays themselves
	int *gather_buffer = NULL;
	if (procno == 0) gather_buffer = (int *) malloc(total_number_of_elements * sizeof(int));
	MPI_Gatherv(
			my_elements, // Items to gather
			my_number_of_elements, // Size of items per process (varied)
			MPI_INT, // Type
			gather_buffer, // Where to put result
			size_buffer, // Size to expect incoming per process (varied)
			displ_buffer, // Offsets of each incoming buffer
			MPI_INT, // Type
			0, // Process number at which to put the result
			comm // Communicator
		   );

	// Peint the gathered array
	if (procno == 0) {
		printf("Gathered:");
		for (int i = 0; i < total_number_of_elements; i++) {
			printf(" %d", gather_buffer[i]);
		}
		printf("\n");
	}

	MPI_Finalize();

	return 0;
}
