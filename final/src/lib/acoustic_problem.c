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
#include <stdlib.h>
#include <stdio.h>
#include <sys/stat.h>

PetscErrorCode solve_case(MPI_Comm comm, PetscMPIInt procno, PetscMPIInt nprocs, struct AcousticCase ac) {
	PetscFunctionBegin;
	// Validate entries if necessary
	//
	//
	//
	
	// File writing variables
	char *dir = (char *) malloc(100);
	char *path = (char *) malloc(100);
	FILE *file;
	
	// Write case directory
	sprintf(dir, "./%s", ac.op.name);
	if (procno == 0) {
		printf("\tCreating directory %s...\n", dir);
		if (mkdir(dir, S_IRWXU | S_IRWXG | S_IRWXO)) {
			printf("\tDirectory creation failed. Does directory \"%s\" already exist?\n", dir);
			return -1;
		}
	}
	
	// Save case info
	if (procno == 0) {
		PetscCall(_write_metadata(ac));
	}
	
	// Problem setup with PETSc
	PetscCall(PetscPrintf(comm, "\tSetting up PETSc problem...\n"));
	Mat A;
	Vec b, zn, znm1, znm2, prov;
	
	// Build matrix
	PetscCall(_build_matrix(comm, procno, nprocs, ac, &A));
	
	if (ac.op.debug) {
		PetscCall(PetscPrintf(comm, "\t\tWriting \"A\" matrix to file...\n"));
		PetscViewer viewer;
		sprintf(path, "%s/A", dir);
		PetscCall(PetscViewerASCIIOpen(PETSC_COMM_WORLD, path, &viewer));
		PetscCall(MatView(A, viewer));
	}
	
	// Build right hand side
	PetscCall(_build_rhs(comm, procno, nprocs, ac, &b));
	
	if (ac.op.debug) {
		PetscCall(PetscPrintf(comm, "\t\tWriting \"b\" vector to file...\n"));
		PetscViewer viewer;
		sprintf(path, "%s/b", dir);
		PetscCall(PetscViewerASCIIOpen(PETSC_COMM_WORLD, path, &viewer));
		PetscCall(VecView(b, viewer));
	}
	
	// Initialize solution variables
	// zn
	PetscCall(_build_z0(comm, procno, nprocs, ac, &zn));
	
	if (ac.op.debug) {
		PetscCall(PetscPrintf(comm, "\t\tWriting \"z0\" vector to file...\n"));
		PetscViewer viewer;
		sprintf(path, "%s/z0", dir);
		PetscCall(PetscViewerASCIIOpen(PETSC_COMM_WORLD, path, &viewer));
		PetscCall(VecView(zn, viewer));
	}
	
	// znm1 and znm2 (empty vectors)
	PetscInt vec_size = ac.sd.nx * ac.sd.ny * 3;
	PetscCall(VecCreate(comm, &znm1));
	PetscCall(VecSetType(znm1, VECSTANDARD)); // Or VECMPI
	PetscCall(VecSetSizes(znm1, PETSC_DECIDE, vec_size));
	PetscCall(VecCreate(comm, &znm2));
	PetscCall(VecSetType(znm2, VECSTANDARD)); // Or VECMPI
	PetscCall(VecSetSizes(znm2, PETSC_DECIDE, vec_size));
	PetscCall(VecCreate(comm, &prov));
	PetscCall(VecSetType(prov, VECSTANDARD)); // Or VECMPI
	PetscCall(VecSetSizes(prov, PETSC_DECIDE, vec_size));
	
	// Starting time loop
	PetscCall(PetscPrintf(comm, "\nStarting time loop...\n"));
	double t = ac.td.ta;
	double dt = (ac.td.tb - ac.td.ta) / ac.td.nt;
	int n = 0;
	PetscCall(PetscPrintf(comm, "\tTime = %8.3f (%6.2f%%)\n", t, n / (double) ac.td.nt * 100));
	// Shuffle solution vectors
	// N/A
	// Calculate latest (except for values at boundaries)
	// N/A
	// Apply boundary conditions to latest
	PetscCall(_apply_dbc(comm, procno, nprocs, ac, zn, t));
	// Write results
	if (n % ac.op.write_every == 0 && (ac.op.write_single || ac.op.write_split)) {
		PetscCall(_write_step(comm, procno, nprocs, ac, zn, n, t));
	}
	if (ac.op.debug) {
		printf("\t\tProcess %d finished with step %d\n", procno, n);
	}
	
	// First order initial step
	n = 1;
	t += dt;
	PetscCall(PetscPrintf(comm, "\tTime = %8.3f (%6.2f%%)\n", t, n / (double) ac.td.nt * 100));
	// Shuffle solution vectors
	PetscCall(VecCopy(zn, znm1));
	// Calculate latest (except for values at boundaries)
	PetscCall(PetscPrintf(comm, "\t\tCalculating interior explicitly...\n"));
	// zn = znm1 + k * ((A * znm1) + b)
	// prov = A * znm1
	PetscCall(MatMult(
		// y = A*x
		A, // A
		znm1, // x
		prov // y
	));
	// prov = prov + b
	PetscCall(VecAXPY(
		// y = a*x + y
		prov, // y
		1, // a
		b // x
	));
	// zn = znm1 + k * prov
	PetscCall(VecAXPBYPCZ(
		// z = a*x + b*y + c*z
		zn, // z
		1, // a
		dt, // b
		0, // c
		znm1, // x
		prov // y
	));
	// Apply boundary conditions to latest
	PetscCall(_apply_dbc(comm, procno, nprocs, ac, zn, t));
	// Write results
	if (n % ac.op.write_every == 0 && (ac.op.write_single || ac.op.write_split)) {
		PetscCall(_write_step(comm, procno, nprocs, ac, zn, n, t));
	}
	if (ac.op.debug) {
		printf("\t\tProcess %d finished with step %d\n", procno, n);
	}
	
	// Second order time loop
	for (n = 2; n <= ac.td.nt; n++) {
		t += dt;
		PetscCall(PetscPrintf(comm, "\tTime = %8.3f (%6.2f%%)\n", t, n / (double) ac.td.nt * 100));
		// Shuffle solution vectors
		PetscCall(VecCopy(znm1, znm2));
		PetscCall(VecCopy(zn, znm1));
		// Calculate latest (except for values at boundaries)
		PetscCall(PetscPrintf(comm, "\t\tCalculating interior explicitly...\n"));
		// zn = znm2 + 2 * k * ((A * znm1) + b)
		// prov = A * znm1
		PetscCall(MatMult(
			// y = A*x
			A, // A
			znm1, // x
			prov // y
		));
		// prov = prov + b
		PetscCall(VecAXPY(
			// y = a*x + y
			prov, // y
			1, // a
			b // x
		));
		// zn = znm2 + 2 * k * prov
		PetscCall(VecAXPBYPCZ(
			// z = a*x + b*y + c*z
			zn, // z
			1, // a
			2 * dt, // b
			0, // c
			znm2, // x
			prov // y
		));
		// Apply boundary conditions to latest
		PetscCall(_apply_dbc(comm, procno, nprocs, ac, zn, t));
		// Write results
		if (n % ac.op.write_every == 0 && (ac.op.write_single || ac.op.write_split)) {
			PetscCall(_write_step(comm, procno, nprocs, ac, zn, n, t));
		}
		if (ac.op.debug) {
			printf("\t\tProcess %d finished with step %d\n", procno, n);
		}
	}

	free(dir);
	free(path);
	
	PetscCall(MatDestroy(&A));
	PetscCall(VecDestroy(&b));
	PetscCall(VecDestroy(&zn));
	PetscCall(VecDestroy(&znm1));
	PetscCall(VecDestroy(&znm2));

	PetscCall(PetscPrintf(comm, "\nCase complete.\n"));
	
	PetscFunctionReturn(0);
}

PetscErrorCode _build_matrix(MPI_Comm comm, PetscMPIInt procno, PetscMPIInt nprocs, struct AcousticCase ac, Mat *rA) {
	// Note here that Mat *rA is a pointer
	PetscFunctionBegin;
	PetscCall(PetscPrintf(comm, "\t\tAssembling \"A\" matrix...\n"));

	// Basic matrix creation
	Mat A;
	PetscCall(MatCreate(comm, &A));
	PetscCall(MatSetType(A, MATAIJ)); // Or MATMPIAIJ
	PetscInt matrix_size = ac.sd.nx * ac.sd.ny * 3;
	PetscCall(MatSetSizes(A, PETSC_DECIDE, PETSC_DECIDE, matrix_size, matrix_size));

	// Loop through domain, bc, and set values by adding
	if (procno == 0) {
		printf("\t\t\tValue assignment loop in serial, for now\n");
		PetscInt i, j, eqn, var;
		PetscScalar v, x, y;
		PetscInt nx = ac.sd.nx;
		PetscInt ny = ac.sd.ny;
		PetscScalar dx = (ac.sd.xb - ac.sd.xa) / (ac.sd.nx - 1);
		PetscScalar dy = (ac.sd.yb - ac.sd.ya) / (ac.sd.ny - 1);
		for (i = 0; i < nx; i++) {
			for (j = 0; j < ny; j++) {
				// SET X AND Y HERE
				x = ac.sd.xa + i * dx;
				y = ac.sd.ya + j * dy;
				// 3 equations associated with DOF (i, j, *)
				PetscInt eij_k0 = (ny * 3) * i + 3 * j; // rho_p
				PetscInt eij_k1 = (ny * 3) * i + 3 * j + 1; // u_p
				PetscInt eij_k2 = (ny * 3) * i + 3 * j + 2; // v_p
				// DOF (i - 1, j, *)
				PetscInt eim1j_k0 = (ny * 3) * (i - 1) + 3 * j; // rho_p
				PetscInt eim1j_k1 = (ny * 3) * (i - 1) + 3 * j + 1; // u_p
				PetscInt eim1j_k2 = (ny * 3) * (i - 1) + 3 * j + 2; // v_p
				// DOF (i + 1, j, *)
				PetscInt eip1j_k0 = (ny * 3) * (i + 1) + 3 * j; // rho_p
				PetscInt eip1j_k1 = (ny * 3) * (i + 1) + 3 * j + 1; // u_p
				PetscInt eip1j_k2 = (ny * 3) * (i + 1) + 3 * j + 2; // v_p
				// DOF (i, j - 1, *)
				PetscInt eijm1_k0 = (ny * 3) * i + 3 * (j - 1); // rho_p
				PetscInt eijm1_k1 = (ny * 3) * i + 3 * (j - 1) + 1; // u_p
				PetscInt eijm1_k2 = (ny * 3) * i + 3 * (j - 1) + 2; // v_p
				// DOF (i, j + 1, *)
				PetscInt eijp1_k0 = (ny * 3) * i + 3 * (j + 1); // rho_p
				PetscInt eijp1_k1 = (ny * 3) * i + 3 * (j + 1) + 1; // u_p
				PetscInt eijp1_k2 = (ny * 3) * i + 3 * (j + 1) + 2; // v_p
			
				// Example
				// eqn = eij_k0;
				// var = eij_k0;
				// v = 1.0;
				// PetscCall(MatSetValues(A, 1, &eqn, 1, &var, &v, ADD_VALUES));
			
				// If on boundary, just Dirichlet condition (zero entry for explicit)
				if (!(i == 0 || i == nx - 1 || j == 0 || j == ny - 1)) {
					/*
					 * Equation for rho_p
					 */
					eqn = eij_k0;
				
					// add -rho_bar * div((u_p, v_p))
					// --> -rho_bar / (2 * dx) * u_p(i + 1, j)
					var = eip1j_k1;
					v = -ac.mp.rho_bar / (2 * dx);
					PetscCall(MatSetValues(A, 1, &eqn, 1, &var, &v, ADD_VALUES));
					// -->  rho_bar / (2 * dx) * u_p(i - 1, j)
					var = eim1j_k1;
					v = ac.mp.rho_bar / (2 * dx);
					PetscCall(MatSetValues(A, 1, &eqn, 1, &var, &v, ADD_VALUES));
					// --> -rho_bar / (2 * dy) * v_p(i, j + 1)
					var = eijp1_k2;
					v = -ac.mp.rho_bar / (2 * dy);
					PetscCall(MatSetValues(A, 1, &eqn, 1, &var, &v, ADD_VALUES));
					// -->  rho_bar / (2 * dy) * v_p(i, j - 1)
					var = eijm1_k2;
					v = ac.mp.rho_bar / (2 * dy);
					PetscCall(MatSetValues(A, 1, &eqn, 1, &var, &v, ADD_VALUES));
				
					// add -dot((u_bar, v_bar), grad(rho_p))
					// --> -u_bar / (2 * dx) * rho_p(i + 1, j)
					var = eip1j_k0;
					v = -ac.mp.u_bar(x, y) / (2 * dx);
					PetscCall(MatSetValues(A, 1, &eqn, 1, &var, &v, ADD_VALUES));
					// -->  u_bar / (2 * dx) * rho_p(i - 1, j)
					var = eim1j_k0;
					v = ac.mp.u_bar(x, y) / (2 * dx);
					PetscCall(MatSetValues(A, 1, &eqn, 1, &var, &v, ADD_VALUES));
					// --> -v_bar / (2 * dy) * rho_p(i, j + 1)
					var = eijp1_k0;
					v = -ac.mp.v_bar(x, y) / (2 * dy);
					PetscCall(MatSetValues(A, 1, &eqn, 1, &var, &v, ADD_VALUES));
					// -->  v_bar / (2 * dy) * rho_p(i, j - 1)
					var = eijm1_k0;
					v = ac.mp.v_bar(x, y) / (2 * dy);
					PetscCall(MatSetValues(A, 1, &eqn, 1, &var, &v, ADD_VALUES));
				
					/*
					 * Equation for u_p
					 */
					eqn = eij_k1;
				
					// add -dot((u_bar, v_bar), (d[u_p]/dx, d[u_p]/dy))
					// >   -u_bar * d[u_p]/dx
					// --> -u_bar / (2 * dx) * u_p(i + 1, j)
					var = eip1j_k1;
					v = -ac.mp.u_bar(x, y) / (2 * dx);
					PetscCall(MatSetValues(A, 1, &eqn, 1, &var, &v, ADD_VALUES));
					// -->  u_bar / (2 * dx) * u_p(i - 1, j)
					var = eim1j_k1;
					v = ac.mp.u_bar(x, y) / (2 * dx);
					PetscCall(MatSetValues(A, 1, &eqn, 1, &var, &v, ADD_VALUES));
					// >   -v_bar * d[u_p]/dy
					// --> -v_bar / (2 * dy) * u_p(i, j + 1)
					var = eijp1_k1;
					v = -ac.mp.v_bar(x, y) / (2 * dy);
					PetscCall(MatSetValues(A, 1, &eqn, 1, &var, &v, ADD_VALUES));
					// -->  v_bar / (2 * dy) * u_p(i, j - 1)
					var = eijm1_k1;
					v = ac.mp.v_bar(x, y) / (2 * dy);
					PetscCall(MatSetValues(A, 1, &eqn, 1, &var, &v, ADD_VALUES));
				
					// add -c^2 / rho_bar * d[rho_p]/dx
					// --> -c * c / (rho_bar * 2 * dx) * rho_p(i + 1, j)
					var = eip1j_k0;
					v = -ac.mp.c * ac.mp.c / (ac.mp.rho_bar * 2 * dx);
					PetscCall(MatSetValues(A, 1, &eqn, 1, &var, &v, ADD_VALUES));
					// --> c * c / (rho_bar * 2 * dx) * rho_p(i - 1, j)
					var = eim1j_k0;
					v = ac.mp.c * ac.mp.c / (ac.mp.rho_bar * 2 * dx);
					PetscCall(MatSetValues(A, 1, &eqn, 1, &var, &v, ADD_VALUES));
				
					/*
					 * Equation for v_p
					 */
					eqn = eij_k2;
				
					// add -dot((u_bar, v_bar), (d[v_p]/dx, d[v_p]/dy))
					// >   -u_bar * d[v_p]/dx
					// --> -u_bar / (2 * dx) * v_p(i + 1, j)
					var = eip1j_k2;
					v = -ac.mp.u_bar(x, y) / (2 * dx);
					PetscCall(MatSetValues(A, 1, &eqn, 1, &var, &v, ADD_VALUES));
					// --> u_bar / (2 * dx) * v_p(i - 1, j)
					var = eim1j_k2;
					v = ac.mp.u_bar(x, y) / (2 * dx);
					PetscCall(MatSetValues(A, 1, &eqn, 1, &var, &v, ADD_VALUES));
					// >   -v_bar * d[v_p]/dy
					// --> -v_bar / (2 * dy) * v_p(i, j + 1)
					var = eijp1_k2;
					v = -ac.mp.v_bar(x, y) / (2 * dy);
					PetscCall(MatSetValues(A, 1, &eqn, 1, &var, &v, ADD_VALUES));
					// --> v_bar / (2 * dy) * v_p(i, j - 1)
					var = eijm1_k2;
					v = ac.mp.v_bar(x, y) / (2 * dy);
					PetscCall(MatSetValues(A, 1, &eqn, 1, &var, &v, ADD_VALUES));
				
					// add -c^2 / rho_bar * d[rho_p]/dy
					// --> -c * c / (rho_bar * 2 * dy) * rho_p(i, j + 1)
					var = eijp1_k0;
					v = -ac.mp.c * ac.mp.c / (ac.mp.rho_bar * 2 * dy);
					PetscCall(MatSetValues(A, 1, &eqn, 1, &var, &v, ADD_VALUES));
					// --> c * c / (rho_bar * 2 * dy) * rho_p(i, j - 1)
					var = eijm1_k0;
					v = ac.mp.c * ac.mp.c / (ac.mp.rho_bar * 2 * dy);
					PetscCall(MatSetValues(A, 1, &eqn, 1, &var, &v, ADD_VALUES));
				}
			}
		}
	}
	
	// Set values, assemble
	PetscCall(MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY));
	PetscCall(MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY));
	
	*rA = A;
	
	PetscFunctionReturn(0);
}

PetscErrorCode _build_rhs(MPI_Comm comm, PetscMPIInt procno, PetscMPIInt nprocs, struct AcousticCase ac, Vec *rb) {
	// Note here that *rb is a pointer
	PetscFunctionBegin;
	PetscCall(PetscPrintf(comm, "\t\tAssembling \"b\" vector...\n"));

	// Basic matrix creation
	Vec b;
	PetscCall(VecCreate(comm, &b));
	PetscCall(VecSetType(b, VECSTANDARD)); // Or VECMPI
	PetscInt vec_size = ac.sd.nx * ac.sd.ny * 3;
	PetscCall(VecSetSizes(b, PETSC_DECIDE, vec_size));

	// Loop through domain, bc, and set values by adding
	if (procno == 0) {
		printf("\t\t\tValue assignment loop in serial, for now\n");
		PetscInt i, j, eqn;
		PetscScalar v, x, y;
		PetscInt nx = ac.sd.nx;
		PetscInt ny = ac.sd.ny;
		PetscScalar dx = (ac.sd.xb - ac.sd.xa) / (ac.sd.nx - 1);
		PetscScalar dy = (ac.sd.yb - ac.sd.ya) / (ac.sd.ny - 1);
		for (i = 0; i < nx; i++) {
			for (j = 0; j < ny; j++) {
				// SET X AND Y HERE
				x = ac.sd.xa + i * dx;
				y = ac.sd.ya + j * dy;
				// 3 equations associated with DOF (i, j, *)
				PetscInt eij_k0 = (ny * 3) * i + 3 * j; // rho_p
				PetscInt eij_k1 = (ny * 3) * i + 3 * j + 1; // u_p
				PetscInt eij_k2 = (ny * 3) * i + 3 * j + 2; // v_p
			
				// Example
				// eqn = eij_k0;
				// v = 1.0;
				// PetscCall(VecSetValue(b, &eqn, &v, ADD_VALUES));
		
				// If on boundary, just Dirichlet condition (zero entry for explicit)
				if (!(i == 0 || i == nx - 1 || j == 0 || j == ny - 1)) {
					/*
					 * Equation for rho_p
					 */
					eqn = eij_k0;
			
					// NADA
		
					/*
					 * Equation for u_p
					 */
					eqn = eij_k1;
			
					// add -dot((u_bar, v_bar), (d[u_bar]/dx, d[u_bar]/dy))
					// >   -u_bar * d[u_bar]/dx
					v = -ac.mp.u_bar(x, y) * (ac.mp.u_bar(x + dx, y) - ac.mp.u_bar(x - dx, y)) / (2 * dx);
					PetscCall(VecSetValue(b, eqn, v, ADD_VALUES));
					// >   -v_bar * d[u_p]/dy
					v = -ac.mp.v_bar(x, y) * (ac.mp.u_bar(x, y + dy) - ac.mp.u_bar(x, y - dy)) / (2 * dy);
					PetscCall(VecSetValue(b, eqn, v, ADD_VALUES));
			
					/*
					 * Equation for v_p
					 */
					eqn = eij_k2;
			
					// add -dot((u_bar, v_bar), (d[v_bar]/dx, d[v_bar]/dy))
					// >   -u_bar * d[v_bar]/dx
					v = -ac.mp.u_bar(x, y) * (ac.mp.v_bar(x + dx, y) - ac.mp.v_bar(x - dx, y)) / (2 * dx);
					PetscCall(VecSetValue(b, eqn, v, ADD_VALUES));
					// >   -v_bar * d[v_bar]/dy
					v = -ac.mp.v_bar(x, y) * (ac.mp.v_bar(x, y + dy) - ac.mp.v_bar(x, y - dy)) / (2 * dy);
					PetscCall(VecSetValue(b, eqn, v, ADD_VALUES));
				}
			}
		}
	}
	
	// Set values, assemble
	PetscCall(VecAssemblyBegin(b));
	PetscCall(VecAssemblyEnd(b));
	
	*rb = b;
	
	PetscFunctionReturn(0);
}

PetscErrorCode _build_z0(MPI_Comm comm, PetscMPIInt procno, PetscMPIInt nprocs, struct AcousticCase ac, Vec *rz) {
	// Note here that *rb is a pointer
	PetscFunctionBegin;
	PetscCall(PetscPrintf(comm, "\t\tAssembling \"z0\" vector...\n"));

	// Basic matrix creation
	Vec z;
	PetscCall(VecCreate(comm, &z));
	PetscCall(VecSetType(z, VECSTANDARD)); // Or VECMPI
	PetscInt vec_size = ac.sd.nx * ac.sd.ny * 3;
	PetscCall(VecSetSizes(z, PETSC_DECIDE, vec_size));

	// Loop through domain, bc, and set values by adding
	if (procno == 0) {
		printf("\t\t\tValue assignment loop in serial, for now\n");
		PetscInt i, j, eqn;
		PetscScalar v, x, y;
		PetscInt nx = ac.sd.nx;
		PetscInt ny = ac.sd.ny;
		PetscScalar dx = (ac.sd.xb - ac.sd.xa) / (ac.sd.nx - 1);
		PetscScalar dy = (ac.sd.yb - ac.sd.ya) / (ac.sd.ny - 1);
		for (i = 0; i < nx; i++) {
			for (j = 0; j < ny; j++) {
				// SET X AND Y HERE
				x = ac.sd.xa + i * dx;
				y = ac.sd.ya + j * dy;
				// 3 equations associated with DOF (i, j, *)
				PetscInt eij_k0 = (ny * 3) * i + 3 * j; // rho_p
				PetscInt eij_k1 = (ny * 3) * i + 3 * j + 1; // u_p
				PetscInt eij_k2 = (ny * 3) * i + 3 * j + 2; // v_p
			
				// Example
				// eqn = eij_k0;
				// v = 1.0;
				// PetscCall(VecSetValue(z, &eqn, &v, ADD_VALUES));
		
				// If on boundary, just Dirichlet condition (zero entry for explicit)
				if (!(i == 0 || i == nx - 1 || j == 0 || j == ny - 1)) {
					/*
					 * Equation for rho_p
					 */
					eqn = eij_k0;
					v = ac.ic.rho_p_0(x, y);
					PetscCall(VecSetValue(z, eqn, v, INSERT_VALUES));
		
					/*
					 * Equation for u_p
					 */
					eqn = eij_k1;
					v = ac.ic.u_p_0(x, y);
					PetscCall(VecSetValue(z, eqn, v, INSERT_VALUES));
			
					/*
					 * Equation for v_p
					 */
					eqn = eij_k2;
					v = ac.ic.v_p_0(x, y);
					PetscCall(VecSetValue(z, eqn, v, INSERT_VALUES));
				}
			}
		}
	}
	
	// Set values, assemble
	PetscCall(VecAssemblyBegin(z));
	PetscCall(VecAssemblyEnd(z));
	
	*rz = z;
	
	PetscFunctionReturn(0);
}

PetscErrorCode _apply_dbc(MPI_Comm comm, PetscMPIInt procno, PetscMPIInt nprocs, struct AcousticCase ac, Vec z, double t) {
	PetscFunctionBegin;
	PetscCall(PetscPrintf(comm, "\t\tApplying Dirichlet boundary conditions...\n"));
	
	// Loop through domain, bc, and set values by adding
	if (procno == 0) {
		printf("\t\t\tValue assignment loop in serial, for now\n");
		PetscInt i, j, eqn;
		PetscScalar v, x, y;
		PetscInt nx = ac.sd.nx;
		PetscInt ny = ac.sd.ny;
		PetscScalar dx = (ac.sd.xb - ac.sd.xa) / (ac.sd.nx - 1);
		PetscScalar dy = (ac.sd.yb - ac.sd.ya) / (ac.sd.ny - 1);
		for (i = 0; i < nx; i++) {
			for (j = 0; j < ny; j++) {
				// SET X AND Y HERE
				x = ac.sd.xa + i * dx;
				y = ac.sd.ya + j * dy;
				// 3 equations associated with DOF (i, j, *)
				PetscInt eij_k0 = (ny * 3) * i + 3 * j; // rho_p
				PetscInt eij_k1 = (ny * 3) * i + 3 * j + 1; // u_p
				PetscInt eij_k2 = (ny * 3) * i + 3 * j + 2; // v_p
			
				// Example
				// eqn = eij_k0;
				// v = 1.0;
				// PetscCall(VecSetValue(z, eqn, v, ADD_VALUES));
		
				// If on boundary, apply Dirichlet condition
				if (i == 0) {
					eqn = eij_k0;
					v = ac.bc.rho_p_xa(x, y, t);
					PetscCall(VecSetValue(z, eqn, v, INSERT_VALUES));

					eqn = eij_k1;
					v = ac.bc.u_p_xa(x, y, t);
					PetscCall(VecSetValue(z, eqn, v, INSERT_VALUES));

					eqn = eij_k2;
					v = ac.bc.v_p_xa(x, y, t);
					PetscCall(VecSetValue(z, eqn, v, INSERT_VALUES));
				} else if (i == nx - 1) {
					eqn = eij_k0;
					v = ac.bc.rho_p_xb(x, y, t);
					PetscCall(VecSetValue(z, eqn, v, INSERT_VALUES));

					eqn = eij_k1;
					v = ac.bc.u_p_xb(x, y, t);
					PetscCall(VecSetValue(z, eqn, v, INSERT_VALUES));

					eqn = eij_k2;
					v = ac.bc.v_p_xb(x, y, t);
					PetscCall(VecSetValue(z, eqn, v, INSERT_VALUES));
				} else if (j == 0) {
					eqn = eij_k0;
					v = ac.bc.rho_p_ya(x, y, t);
					PetscCall(VecSetValue(z, eqn, v, INSERT_VALUES));

					eqn = eij_k1;
					v = ac.bc.u_p_ya(x, y, t);
					PetscCall(VecSetValue(z, eqn, v, INSERT_VALUES));

					eqn = eij_k2;
					v = ac.bc.v_p_ya(x, y, t);
					PetscCall(VecSetValue(z, eqn, v, INSERT_VALUES));
				} else if (j == ny - 1) {
					eqn = eij_k0;
					v = ac.bc.rho_p_yb(x, y, t);
					PetscCall(VecSetValue(z, eqn, v, INSERT_VALUES));

					eqn = eij_k1;
					v = ac.bc.u_p_yb(x, y, t);
					PetscCall(VecSetValue(z, eqn, v, INSERT_VALUES));

					eqn = eij_k2;
					v = ac.bc.v_p_yb(x, y, t);
					PetscCall(VecSetValue(z, eqn, v, INSERT_VALUES));
				}
			}
		}
	}
	
	// Set values, assemble
	PetscCall(VecAssemblyBegin(z));
	PetscCall(VecAssemblyEnd(z));
	
	PetscFunctionReturn(0);
}

PetscErrorCode _write_metadata(struct AcousticCase ac) {
	PetscFunctionBegin;
	printf("\tSaving case metadata...\n");
	
	char *dir = (char *) malloc(100);
	sprintf(dir, "./%s", ac.op.name);
	
	char *path = (char *) malloc(100);
	FILE *file;
	
	sprintf(path, "%s/metadata.json", dir);
	file = fopen(path, "w");
	
	// SpaceDomain
	fprintf(file, "{\n\t\"SpaceDomain\": {\n");
	fprintf(file, "\t\t\"dim\": 2,\n");
	fprintf(file, "\t\t\"xa\": %f,\n", ac.sd.xa);
	fprintf(file, "\t\t\"xb\": %f,\n", ac.sd.xb);
	fprintf(file, "\t\t\"nx\": %d,\n", ac.sd.nx);
	fprintf(file, "\t\t\"ya\": %f,\n", ac.sd.ya);
	fprintf(file, "\t\t\"yb\": %f,\n", ac.sd.yb);
	fprintf(file, "\t\t\"ny\": %d\n", ac.sd.ny);
	fprintf(file, "\t},\n");
	
	// TimeDomain
	fprintf(file, "\t\"TimeDomain\": {\n");
	fprintf(file, "\t\t\"dim\": 2,\n");
	fprintf(file, "\t\t\"xa\": %f,\n", ac.sd.xa);
	fprintf(file, "\t\t\"xb\": %f,\n", ac.sd.xb);
	fprintf(file, "\t\t\"nx\": %d,\n", ac.sd.nx);
	fprintf(file, "\t\t\"ya\": %f,\n", ac.sd.ya);
	fprintf(file, "\t\t\"yb\": %f,\n", ac.sd.yb);
	fprintf(file, "\t\t\"ny\": %d\n", ac.sd.ny);
	fprintf(file, "\t},\n");
	
	// OutputParameters
	fprintf(file, "\t\"OutputParameters\": {\n");
	fprintf(file, "\t\t\"name\": \"%s\",\n", ac.op.name);
	fprintf(file, "\t\t\"write_every\": %d,\n", ac.op.write_every);
	if (ac.op.debug) {
		fprintf(file, "\t\t\"debug\": true,\n");
	} else {
		fprintf(file, "\t\t\"debug\": false,\n");
	}
	if (ac.op.write_single) {
		fprintf(file, "\t\t\"write_single\": true,\n");
	} else {
		fprintf(file, "\t\t\"write_single\": false,\n");
	}
	if (ac.op.write_split) {
		fprintf(file, "\t\t\"write_split\": true\n");
	} else {
		fprintf(file, "\t\t\"write_split\": false\n");
	}
	fprintf(file, "\t}\n}\n");
	
	fclose(file);
	
	PetscFunctionReturn(0);
}

PetscErrorCode _write_step(MPI_Comm comm, PetscMPIInt procno, PetscMPIInt nprocs, struct AcousticCase ac, Vec z, int n, double t) {
	PetscFunctionBegin;
	PetscCall(PetscPrintf(comm, "\t\tWriting to file...\n"));
	
	// Write step subdirectory
	char *dir = (char *) malloc(100);
	sprintf(dir, "./%s/%d", ac.op.name, n);
	if (procno == 0) {
		if (mkdir(dir, S_IRWXU | S_IRWXG | S_IRWXO)) {
			printf("\t\t\tDirectory creation failed. Does directory \"%s\" already exist?\n", dir);
			PetscFunctionReturn(1);
		}
	}
	
	char *path = (char *) malloc(100);
	FILE *file;
	
	// t
	sprintf(path, "%s/t", dir);
	if (procno == 0) {
		file = fopen(path, "w");
		fprintf(file, "%e\n", t);
		fclose(file);
	}
	
	if (ac.op.write_single) {
		PetscViewer viewer;
		sprintf(path, "%s/z", dir);
		PetscCall(PetscViewerASCIIOpen(PETSC_COMM_WORLD, path, &viewer));
		PetscCall(VecView(z, viewer));
	}
	
	if (procno == 0 && ac.op.write_split) {
		printf("WARNING: Splitting solution vector writing deferred, for now\n");
		/*
		sprintf(path, "%s/rho_p", dir);
		FILE *file_rho_p = fopen(path, "w");
		
		sprintf(path, "%s/u_p", dir);
		FILE *file_u_p = fopen(path, "w");
		
		sprintf(path, "%s/v_p", dir);
		FILE *file_v_p = fopen(path, "w");
		
		PetscScalar v;
		for (int i = 0; i < ac.sd.nx; i++) {
			for (int j = 0; j < ac.sd.ny; j++) {
				// 3 equations associated with DOF (i, j, *)
				PetscInt eij_k0 = (ac.sd.ny * 3) * i + 3 * j; // rho_p
				PetscInt eij_k1 = (ac.sd.ny * 3) * i + 3 * j + 1; // u_p
				PetscInt eij_k2 = (ac.sd.ny * 3) * i + 3 * j + 2; // v_p
				
				PetscCall(VecGetValues(z, 1, &eij_k0, &v));
				fprintf(file_rho_p, "%e\n", v);
				
				PetscCall(VecGetValues(z, 1, &eij_k1, &v));
				fprintf(file_u_p, "%e\n", v);
				
				PetscCall(VecGetValues(z, 1, &eij_k2, &v));
				fprintf(file_v_p, "%e\n", v);
			}
		}
		
		fclose(file_rho_p);
		fclose(file_u_p);
		fclose(file_v_p);
		*/
	}
	
	free(dir);
	free(path);
	
	PetscFunctionReturn(0);
}