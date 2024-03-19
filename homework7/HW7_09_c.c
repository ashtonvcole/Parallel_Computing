#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <omp.h>

int main(int argc, char **argv) {
	int debug = 0;

	long int vectorsize = argc == 2 ? atoi(argv[1]) : 100;
	printf("Using vector size=%ld\n", vectorsize);
	
	double diag = 2.01;
	double *solution = (double *) malloc(vectorsize * sizeof(double));
	double *rhs = (double *) malloc(vectorsize * sizeof(double));
	double *xvector = (double *) malloc(vectorsize * sizeof(double));
	double *tvector = (double *) malloc(vectorsize * sizeof(double));

	if (!solution || !rhs || !xvector || !tvector) {
		printf("Allocation failed\n");
		return 1;
	}

	// A(i, i) = diag
	// A(i, i - 1) = A(i, i + 1) = -1
	// Let the solution be a vector of ones
	// Since A is tri-diagoinal in the form above,
	// RHS must be diag - 2
#pragma omp parallel for
	{
		for (long int i = 0; i < vectorsize; i++) {
			solution[i] = 1.;
			rhs[i] = diag - 2.;
			xvector[i] = i * 1. / vectorsize; // Initial guess
							  // Not necessarily good
		}
	}
	rhs[0] = diag - 1.;
	rhs[vectorsize - 1] = diag - 1.;

	double error0, error; // Initial and present RSS error
	double tstart = omp_get_wtime();
	int iteration;

	for(iteration = 0; /* NADA */; iteration++) {
		// Error by root of the sum of squares
		error = 0.;
		
		// Each calculation is independent
		// But the results are all added to the same variable
		// error += ... is the same as
		// error = error + ...
		// There is a potential race condition here
		// Thread a accesses
		// Thread b accesses
		// Thread a updates
		// Thread b updates
#pragma omp parallel for reduction(+:error)
		for (long int i = 0; i < vectorsize; i++) {
			error += pow(xvector[i] - solution[i], 2);
		}
		error = sqrt(error);
		if (debug) printf("[%d] error=%e\n", iteration, error);

		// Exit the loop if error decreases by
		// at least 3 orders of magnitude
		if (iteration == 0) error0 = error;
		else if (error < error0 * 1.e-3) break;

		// Compute x for next iteration
		// Store in a temporary vector "tvector"
		// x[i, n + 1] =
		// A[i, i] * (b[i] - sum(a[i, j] * x[j], j != i))
		// Each iteration is independent
		// nxb is locally defined
#pragma omp parallel for
		for (long int i = 0; i < vectorsize; i++) {
			// Start with b[i]
			double nxb = rhs[i];
			// Subtract -1 * x[*] for off diagonals
			if (i < vectorsize - 1) nxb += xvector[i + 1];
			if (i > 0) nxb += xvector[i - 1];
			// Not sure why division here, but it works...
			tvector[i] = nxb / diag;
		}
		// Save to x
#pragma omp parallel for
		for (long int i = 0; i < vectorsize; i++) {
			xvector[i] = tvector[i];
		}
		
		if (debug) {
			for (long int i = 0; i < vectorsize; i++) {
				printf("%5.2f\n ", xvector[i]);
			}
		}
	}

	// Get elapsed time
	double duration = omp_get_wtime() - tstart;
	printf("Converged in %d iterations to %e reduction in time=%7.3f\n", iteration, error / error0, duration);

	return 0;
}
