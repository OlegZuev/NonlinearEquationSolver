#pragma once
#include <iostream>
#define PRECISION 6

const int variant = 7;
const int n = 2;
const double eps = 1E-4;

class Vector;
void print_header(std::ostream& ostr);

double error_estimate(double q, double x, double next_x, double y, double next_y);

void set_x0_y0(double& x, double& y);

double func1(double x, double y);

double func2(double x, double y);

double compute_next_x(double y);

double compute_next_y(double x);

void simple_iteration_method(double x, double y, std::ostream& ostr);

void newton_method(double x, double y, std::ostream& ostr);

void gradient_descent_method(double x, double y, std::ostream& ostr);

void print_computed_jacobian(std::ostream& ostr);

void print_nonlinear_system(std::ostream& ostr);

void print_jacobian(std::ostream& ostr);

void print_derivatives_matrix(std::ostream& ostr);

double func1_derivative_x(double x, double y);

double func1_derivative_y(double x, double y);

double func2_derivative_x(double x, double y);

double func2_derivative_y(double x, double y);