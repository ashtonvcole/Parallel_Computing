/*
 * HW6.c
 * Stride-scatter
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

	if (nprocs < 2) {
		printf("This program needs at least two processes\n");
	}

	int sender = 0;
	int localsize = 10;

	int *data = NULL; // Only on root

	MPI_Datatype scattertype, interleavetype;

	if (procno == sender) {
		// Create data array
		int ndata = localsize * nprocs;
		data = (int *) malloc(ndata * sizeof(int));
		if (!data) {
			printf("Out of memory\n");
			MPI_Abort(comm, 0);
		}
		for (int i = 0; i < ndata; i++) {
			data[i] = i;
		}

		// Create strided data type
		int count, stride, blocklength;
		MPI_Aint l, e;
		count = localsize;
		stride = nprocs;
		blocklength = 1;
		MPI_Type_vector(count, blocklength, stride, MPI_INT, &scattertype);
		MPI_Type_commit(&scattertype);
		
		// Derive interleavedtype from scatterype
		MPI_Type_get_extent(scattertype, &l, &e);
		e = blocklength * sizeof(int);
		MPI_Type_create_resized(scattertype, l, e, &interleavetype);
		MPI_Type_commit(&interleavetype);
	}

	// Send out data
	int *mydata = (int *) malloc(localsize * sizeof(int));
	MPI_Scatter(
			data,
			1,
			interleavetype,
			mydata,
			localsize,
			MPI_INT,
			sender,
			comm
		   );

	// Verify
	for (int i = 0; i < localsize; i++) {
		if (mydata[i] % nprocs != procno) {
			printf("[%d] received element=%d, should be %d\n", procno, mydata[i], i * nprocs + procno);
		}
	}

	if (procno == 0) printf("Finished\n");

	if (procno == sender) {
		MPI_Type_free(&interleavetype);
		MPI_Type_free(&scattertype);
		free(data);
	}

	MPI_Finalize();

	return 0;
}
