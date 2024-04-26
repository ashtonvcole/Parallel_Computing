/*
 * Linear Acoustic Perturbation Solver
 * Final Project for Parallel Computing Class
 *
 * "acoustic_problem.h"
 * A header file defining public functions and structures for the solver.
 * 
 * (C) 2024 Ashton Cole. All rights reserved.
 */

struct SpaceDomain {
	double xa; // Left domain boundary
	double xb; // Right domain boundary
	int nx; // Number of solution nodes, boundary-inclusive
	double ya; // Bottom domain boundary
	double yb; // Top domain boundary
	int ny; // Number of solution nodes, boundary-inclusive
};

struct TimeDomain {
	double ta; // Start time
	double tb; // End time
	int nt; // Number of time steps calculated (excludes start)
};

struct MediumProperties {
	double c; // Speed of sound
	double rho_bar; // Density
	double (*u_bar)(double x, double y); // Horizontal-component velocity profile
	double (*v_bar)(double x, double y); // Vertical-component velocity profile
};

struct InitialConditions {
	double (*rho_p_0)(double x, double y); // Acoustic perturbation
	double (*u_p_0)(double x, double y); // Horizontal-component velocity perturbation
	double (*v_p_0)(double x, double y); // Vertical-component velocity perturbation
};

struct BoundaryConditions {
	double (*rho_p_xa)(double x, double y, double t);
	double (*rho_p_xb)(double x, double y, double t);
	double (*rho_p_ya)(double x, double y, double t);
	double (*rho_p_yb)(double x, double y, double t);
	double (*u_p_xa)(double x, double y, double t);
	double (*u_p_xb)(double x, double y, double t);
	double (*u_p_ya)(double x, double y, double t);
	double (*u_p_yb)(double x, double y, double t);
	double (*v_p_xa)(double x, double y, double t);
	double (*v_p_xb)(double x, double y, double t);
	double (*v_p_ya)(double x, double y, double t);
	double (*v_p_yb)(double x, double y, double t);
};

struct OutputParameters {
	char *name; // Directory name
	int write_every; // Write every _th time step
};

struct AcousticCase {
	struct SpaceDomain sd;
	struct TimeDomain td;
	struct MediumProperties mp;
	struct InitialConditions ic;
	struct BoundaryConditions bc;
	struct OutputParameters op;
	// add maps ijk to e and vice-versa?
};

int solve_case(struct AcousticCase ac);
int _build_matrix();
int _build_rhs();
int _solve_matrix_system();
int _write_metadata(struct AcousticCase ac);
int _write_step(struct AcousticCase ac, int n);