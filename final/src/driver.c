/*
 * Linear Acoustic Perturbation Solver
 * Final Project for Parallel Computing Class
 *
 * "driver.c"
 * A main file used to generate an executable.
 * 
 * (C) 2024 Ashton Cole. All rights reserved.
 */

#include "acoustic_problem.h"
#include "driver_functions.h"
#include <stdlib.h>
#include <stdio.h>

int main(int argc, char **argv) {
	printf("Linear Acoustic Perturbation Solver\n");
	printf("Final Project for Parallel Computing Class\n");
	printf("(C) 2024 Ashton Cole. All rights reserved.\n");
	printf("\nStarting driver script...\n");
	
	struct SpaceDomain sd;
	sd.xa = 0; // Left domain boundary
	sd.xb = 10; // Right domain boundary
	sd.nx = 1001; // Number of solution nodes, boundary-inclusive
	sd.ya = 0; // Bottom domain boundary
	sd.yb = 1; // Top domain boundary
	sd.ny = 101; // Number of solution nodes, boundary-inclusive
	
	struct TimeDomain td;
	td.ta = 0; // Start time
	td.tb = 10; // End time
	td.nt = 1000; // Number of time steps calculated (excludes start)
	
	struct MediumProperties mp;
	mp.c = 1; // Speed of sound
	mp.rho_bar = 1; // Density
	mp.u_bar = u_bar; // Horizontal-component velocity profile
	mp.v_bar = v_bar; // Vertical-component velocity profile
	
	struct InitialConditions ic;
	ic.rho_p_0 = rho_p_0; // Acoustic perturbation
	ic.u_p_0 = u_p_0; // Horizontal-component velocity perturbation
	ic.v_p_0 = v_p_0; // Vertical-component velocity perturbation
	
	struct BoundaryConditions bc;
	bc.rho_p_xa = rho_p_xa;
	bc.rho_p_xb = rho_p_xb;
	bc.rho_p_ya = rho_p_ya;
	bc.rho_p_yb = rho_p_yb;
	bc.u_p_xa = u_p_xa;
	bc.u_p_xb = u_p_xb;
	bc.u_p_ya = u_p_ya;
	bc.u_p_yb = u_p_yb;
	bc.v_p_xa = v_p_xa;
	bc.v_p_xb = v_p_xb;
	bc.v_p_ya = v_p_ya;
	bc.v_p_yb = v_p_yb;
	
	struct OutputParameters op;
	op.name = "test";
	op.write_every = 2;
	
	struct AcousticCase ac;
	ac.sd = sd;
	ac.td = td;
	ac.mp = mp;
	ac.ic = ic;
	ac.bc = bc;
	ac.op = op;
	
	solve_case(ac);
	
	return 0;
}
