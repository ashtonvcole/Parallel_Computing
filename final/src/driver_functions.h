/*
 * Linear Acoustic Perturbation Solver
 * Final Project for Parallel Computing Class
 *
 * "acoustic_problem.c"
 * A file defining functions used by the solver for a particular case.
 * 
 * (C) 2024 Ashton Cole. All rights reserved.
 */

#define PI 3.14159265358979323846

// Velocity profile
double u_bar(double x, double y);
double v_bar(double x, double y);

// Initial conditions
double rho_p_0(double x, double y);
double u_p_0(double x, double y);
double v_p_0(double x, double y);

// Boundary conditions
double rho_p_xa(double x, double y, double t);
double rho_p_xb(double x, double y, double t);
double rho_p_ya(double x, double y, double t);
double rho_p_yb(double x, double y, double t);
double u_p_xa(double x, double y, double t);
double u_p_xb(double x, double y, double t);
double u_p_ya(double x, double y, double t);
double u_p_yb(double x, double y, double t);
double v_p_xa(double x, double y, double t);
double v_p_xb(double x, double y, double t);
double v_p_ya(double x, double y, double t);
double v_p_yb(double x, double y, double t);