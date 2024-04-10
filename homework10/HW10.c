/*
 * HW10.c
 * Vector Sum
 * Ashton Cole
 * COE 379L: Parallel Computing
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <omp.h>

int main(int argc, char **argv) {
	int vectorsize = 10;
	if (argc == 2) vectorsize = atoi(argv[1]);
	printf("Using vectorsize: %d\n", vectorsize);

	int nthreads;

#pragma omp parallel
#pragma omp master
	nthreads = omp_get_num_threads();

	double *invec = (double *) malloc(vectorsize * sizeof(double));
	double *outvec = (double *) malloc (vectorsize * sizeof(double));

	if (!invec || !outvec) {
		printf("Could not allocate\n");
		return 1;
	}

	for (int i = 0; i < vectorsize; i++) {
		invec[i] = rand() / (double) RAND_MAX;
		outvec[i] = 0;
	}

	const int nloops = 500;
	double *loopcoeff = (double *) malloc(nloops * sizeof(double));
	loopcoeff[0] = 1.;
	for (int iloop = 1; iloop < nloops; iloop++) {
		loopcoeff[iloop] = loopcoeff[iloop - 1] * 1.1;
	}

	{
		double tstart = omp_get_wtime();
		double factor = 1;
		for (int iloop = 0; iloop < nloops; iloop++) {
			for (int i = 0; i < vectorsize; i++) {
				outvec[i] += invec[i] * loopcoeff[iloop];
			}
		}
		double duration = omp_get_wtime() - tstart;
		printf("Sequential t= %8.5f sec\n", duration);
	}

	{
		double tstart = omp_get_wtime();
		for (int iloop = 0; iloop < nloops; iloop++) {
#pragma omp parallel for
			for (int i = 0; i < vectorsize; i++) {
				outvec[i] += invec[i] * loopcoeff[iloop];
			}
		}
		double duration = omp_get_wtime() - tstart;
		printf("Threads %2d t= %8.5f sec\n", nthreads, duration);
	}

	double s = 0.;
	for (int i = 0; i < vectorsize; i++) {
		s += outvec[i];
	}
	if (s < 0) printf("%e\n", s);

	free(invec);
	free(outvec);
	free(loopcoeff);

	return 0;
}
