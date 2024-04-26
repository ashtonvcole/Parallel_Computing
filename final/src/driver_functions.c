/*
 * Linear Acoustic Perturbation Solver
 * Final Project for Parallel Computing Class
 *
 * "acoustic_problem.c"
 * A file defining functions used by the solver for a particular case.
 * 
 * (C) 2024 Ashton Cole. All rights reserved.
 */

#include "driver_functions.h"
#include <math.h>



/*
 * Velocity Profile
 */

double u_bar(double x, double y) {
	return 0.1;
}

double v_bar(double x, double y) {
	return 0;
}



/*
 * Initial Conditions
 */

double rho_p_0(double x, double y) {
	if (1 < x && x < 2) {
		return 0.0001 * 0.5 * sin(2 * PI * (x - 1));
	} else {
		return 0;
	}
}

double u_p_0(double x, double y) {
	if (1 < x && x < 2) {
		return 0.0001 * 0.5 * sin(2 * PI * (x - 1));
	} else {
		return 0;
	}
}

double v_p_0(double x, double y) {
	return 0;
}



/*
 * Boundary Conditions
 */

double rho_p_xa(double x, double y, double t) {
	return 0;
}

double rho_p_xb(double x, double y, double t) {
	return 0;
}

double rho_p_ya(double x, double y, double t) {
	return 0;
}

double rho_p_yb(double x, double y, double t) {
	return 0;
}

double u_p_xa(double x, double y, double t) {
	return 0;
}

double u_p_xb(double x, double y, double t) {
	return 0;
}

double u_p_ya(double x, double y, double t) {
	return 0;
}

double u_p_yb(double x, double y, double t) {
	return 0;
}

double v_p_xa(double x, double y, double t) {
	return 0;
}

double v_p_xb(double x, double y, double t) {
	return 0;
}

double v_p_ya(double x, double y, double t) {
	return 0;
}

double v_p_yb(double x, double y, double t) {
	return 0;
}