/*
 * Linear Acoustic Perturbation Solver
 * Final Project for Parallel Computing Class
 *
 * "acoustic_problem.c"
 * A file defining functions and structures for the solver.
 * 
 * (C) 2024 Ashton Cole. All rights reserved.
 */

#include "acoustic_problem.h"
#include "petscksp.h"
#include <stdlib.h>
#include <stdio.h>
#include <sys/stat.h>

int solve_case(struct AcousticCase ac) {
	printf("\nStarting case solution...\n");
	
	// Validate entries if necessary
	
	// Write case directory
	char *dir = (char *) malloc(100);
	sprintf(dir, "./%s", ac.op.name);
	printf("\tCreating directory %s\n", dir);
	if (mkdir(dir, S_IRWXU | S_IRWXG | S_IRWXO)) {
		printf("\tDirectory creation failed. Does directory \"%s\" already exist?\n", dir);
		return -1;
	}
	
	// Save case info
	
	// Initialize variables
	// - vectors, matrix, bcs
	double t = ac.td.ta;
	double dt = (ac.td.tb - ac.td.ta) / ac.td.nt;
	
	// Time loop
	printf("\nStarting time loop...\n");
	for (int n = 1; n <= ac.td.nt; n++) {
		t += dt;
		printf("\tTime = %8.3f (%6.2f%%)\n", t, n / (double) ac.td.nt * 100);
		// Update system boundary conditions
		// Solve system
		// Write to file based on write_every
		// STILL NEED TO ACTUALLY WRITE DATA
		// MAKE SURE ONLY ON MAIN THREAD
		if (n % ac.op.write_every == 0) {
			if (_write_step(ac, n)) {
				printf("\tFailure in _write_step(), terminating case.\n");
				return -1;
			}
		}
	}
	
	free(dir);
	
	printf("\nCase complete.\n");
	
	return 0;
}

int _build_matrix() {
	return 0;
}

int _build_rhs() {
	return 0;
}

int _solve_matrix_system() {
	return 0;
}

int _write_metadata(struct AcousticCase ac) {
	return 0;
}

int _write_step(struct AcousticCase ac, int n) {
	// Write case directory
	char *dir = (char *) malloc(100);
	sprintf(dir, "./%s/%d", ac.op.name, n);
	printf("\t\tWriting to file\n");
	if (mkdir(dir, S_IRWXU | S_IRWXG | S_IRWXO)) {
		printf("\tDirectory creation failed. Does directory \"%s\" already exist?\n", dir);
		return -1;
	}
	
	FILE *file;
	char *path = malloc(100);
	
	// rho_p
	sprintf(path, "%s/rho_p", dir);
	file = fopen(path, "w");
	fprintf(file, "FILLER DATA");
	fclose(file);
	
	// u_p
	sprintf(path, "%s/u_p", dir);
	file = fopen(path, "w");
	fprintf(file, "FILLER DATA");
	fclose(file);
	
	// v_p
	sprintf(path, "%s/v_p", dir);
	file = fopen(path, "w");
	fprintf(file, "FILLER DATA");
	fclose(file);
	
	free(dir);
	free(path);
	
	return 0;
}