#include "Functions.h"
#include <utility>
#include <iomanip>
#include <cmath>

double jacobian[2][2];
double derivatives_matrix[2][2];

/**
 * Find error estimate for arguments
 */
double error_estimate_third_norm(double x, double next_x, double y, double next_y) {
	double dif_x = abs(x - next_x);
	double dif_y = abs(y - next_y);
	return sqrt(dif_x * dif_x + dif_y * dif_y);
}

/**
 * Compute first norm of jacobian
 */
double get_norm_jacobian() {
	double result = -1;
	for (int i = 0; i < n; ++i) {
		double sum = 0;
		for (int j = 0; j < n; ++j) {
			sum += fabs(jacobian[i][j]);
		}

		if (sum > result) {
			result = sum;
		}
	}

	return result;
}

/**
 * Get into x and y initial value
 */
void set_x0_y0(double& x, double& y) {
	switch (variant) {
	case 5:
		x = 0.5;
		y = -0.1;
		return;
	case 7:
		x = 1.7;
		y = 0.5;
		return;
	case 0:
		x = -0.1;
		y = 0.5;
		return;
	}

	throw std::invalid_argument("Invalid variant");
}

/**
 *	�ompute jacobian in accordance with the function
 */
void compute_jacobian(double x, double y) {
	switch (variant) {
	case 5:
		jacobian[0][0] = cos(x + 0.5);
		jacobian[0][1] = 0;
		jacobian[1][0] = 0;
		jacobian[1][1] = sin(y - 2);
		return;
	case 7:
		jacobian[0][0] = -cos(x - 1);
		jacobian[0][1] = 0;
		jacobian[1][0] = 0;
		jacobian[1][1] = cos(y + 1);
		return;
	case 0:
		jacobian[0][0] = 0;
		jacobian[0][1] = cos(y + 0.5);
		jacobian[1][0] = sin(x - 2);
		jacobian[1][1] = 0;
		return;
	}

	throw std::invalid_argument("Invalid variant");
}

/**
 * �ompute derivatives matrix in accordance with the function
 */
void compute_derivatives_matrix(double x, double y) {
	derivatives_matrix[0][0] = func1_derivative_x(x, y);
	derivatives_matrix[0][1] = func1_derivative_y(x, y);
	derivatives_matrix[1][0] = func2_derivative_x(x, y);
	derivatives_matrix[1][1] = func2_derivative_y(x, y);
}

/**
 * Reverse derivatives matrix
 */
void reverse_derivatives_matrix() {
	double det = derivatives_matrix[0][0] * derivatives_matrix[1][1] - derivatives_matrix[0][1] * derivatives_matrix[1][
		0];
	std::swap(derivatives_matrix[0][0], derivatives_matrix[1][1]);
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			derivatives_matrix[i][j] /= (i + j == 1) ? -det : det;
		}
	}
}

/**
 * First function in system (depending on the variant)
 */
double func1(double x, double y) {
	switch (variant) {
	case 5:
		return sin(x + 0.5) - y - 1; // y = sin(x + 0.5) - 1	
	case 7:
		return sin(x - 1) - 1.3 + y; // y = -sin(x - 1) + 1.3
	case 0:
		return  +1 - sin(y + 0.5) + x;
	}

	throw std::invalid_argument("Invalid variant");
}

/**
 * Derivative of the first function in system by variable x (depending on the variant)
 */
double func1_derivative_x(double x, double y) {
	switch (variant) {
	case 5:
		return cos(x + 0.5); // sin(x + 0.5) - y - 1;	
	case 7:
		return cos(x - 1); //sin(x - 1) - 1.3 + y;
	case 0:
		return 1; //  +1 - sin(y + 0.5) + x;
	}

	throw std::invalid_argument("Invalid variant");
}

/**
 * Derivative of the first function in system by variable y (depending on the variant)
 */
double func1_derivative_y(double x, double y) {
	switch (variant) {
	case 5:
		return -1; // sin(x + 0.5) - y - 1;	
	case 7:
		return 1; //sin(x - 1) - 1.3 + y;
	case 0:
		return -cos(y + 0.5); //   +1 - sin(y + 0.5) + x;
	}

	throw std::invalid_argument("Invalid variant");
}

/**
 * Second function in system (depending on the variant)
 */
double func2(double x, double y) {
	switch (variant) {
	case 5:
		return cos(y - 2) + x; // x = -cos(y - 2)
	case 7:
		return x - sin(y + 1) - 0.8; // x = sin(y + 1) + 0.8
	case 0:
		return cos(x - 2) + y;
	}

	throw std::invalid_argument("Invalid variant");
}

/**
 * Derivative of the second function in system by variable x (depending on the variant)
 */
double func2_derivative_x(double x, double y) {
	switch (variant) {
	case 5:
		return 1; // cos(y - 2) + x;
	case 7:
		return 1; //x - sin(y + 1) - 0.8;
	case 0:
		return -sin(x - 2); //  cos(x - 2) + y;
	}

	throw std::invalid_argument("Invalid variant");
}

/**
 * Derivative of the second function in system by variable y (depending on the variant)
 */
double func2_derivative_y(double x, double y) {
	switch (variant) {
	case 5:
		return -sin(y - 2); // cos(y - 2) + x;
	case 7:
		return -cos(y + 1); //x - sin(y + 1) - 0.8;
	case 0:
		return 1; //  cos(x - 2) + y;
	}

	throw std::invalid_argument("Invalid variant");
}

/**
 * Compute next x for simple iteration method
 */
double compute_next_x(double y) {
	switch (variant) {
	case 5:
		return -cos(y - 2);
	case 7:
		return sin(y + 1) + 0.8;
	case 0:
		return -1 + sin(y + 0.5);
	}

	throw std::invalid_argument("Invalid variant");
}

/**
 * Compute next y for simple iteration method
 */
double compute_next_y(double x) {
	switch (variant) {
	case 5:
		return sin(x + 0.5) - 1;
	case 7:
		return -sin(x - 1) + 1.3;
	case 0:
		return -cos(x - 2);
	}

	throw std::invalid_argument("Invalid variant");
}

/**
 * This function is used in gradient_descent_method for minimizations
 */
double func_for_minimize(double x, double y) {
	return pow(func1(x, y), 2) + pow(func2(x, y), 2);
}

/**
 * Simple iteration method for solving non linear equation
 */
void simple_iteration_method(double x, double y, std::ostream& ostr) {
	double error;
	int itr = 0;

	ostr << "Simple iteration method" << std::endl;
	print_nonlinear_system(ostr);
	ostr << std::endl << "Jacobian:" << std::endl;
	print_jacobian(ostr);

	ostr << std::endl << "Values:" << std::endl;
	compute_jacobian(x, y);
	print_computed_jacobian(ostr);

	//Header
	ostr << "                               Error" << std::endl;
	ostr << "  Itr |    x     |     y     | Estimate     |      F1    |      F2    | Jacobian norm " << std::endl;

	do {
		itr++;

		double next_x = compute_next_x(y);
		double next_y = compute_next_y(x);

		compute_jacobian(next_x, next_y);
		double q = get_norm_jacobian();
		error = error_estimate_third_norm(x, next_x, y, next_y);

		ostr << std::fixed << std::setw(4) << itr << "  | " << std::setprecision(PRECISION) << std::setw(8) << next_x << " | " << std::setw(9) << next_y <<
			" | " << std::setw(10) << std::scientific << sqrt(func1(next_x, next_y) * func1(next_x, next_y) + func2(next_x, next_y) * func2(next_x, next_y)) << " | " << std::setprecision(PRECISION / 2) << std::
			setw(10) << func1(next_x, next_y) << " | " << std::setw(10) << func2(next_x, next_y) << " | " << std::fixed << std::setprecision(PRECISION) << q << std::endl;

		x = next_x;
		y = next_y;
	} while (error > eps);
}

/**
 * Newton method for solving non linear equation
 */
void newton_method(double x, double y, std::ostream& ostr) {
	double error;
	int itr = 0;

	//Header
	ostr << "Newton method" << std::endl << std::endl;
	ostr << "Derivatives matrix:" << std::endl;
	print_derivatives_matrix(ostr);
	ostr << "                               Error" << std::endl;
	ostr << "  Itr |    x     |     y     | Estimate     |      F1    |      F2    " << std::endl;

	do {
		itr++;

		compute_derivatives_matrix(x, y);
		reverse_derivatives_matrix();

		double computed_func1 = func1(x, y);
		double computed_func2 = func2(x, y);

		double next_x = x - (derivatives_matrix[0][0] * computed_func1 + derivatives_matrix[0][1] * computed_func2);
		double next_y = y - (derivatives_matrix[1][0] * computed_func1 + derivatives_matrix[1][1] * computed_func2);

		error = error_estimate_third_norm(x, next_x, y, next_y);

		ostr << std::fixed << std::setw(4) << itr << "  | " << std::setprecision(PRECISION) << std::setw(8) << next_x << " | " << std::setw(9) << next_y <<
			" | " << std::setw(10) << std::scientific << sqrt(func1(next_x, next_y) * func1(next_x, next_y) + func2(next_x, next_y) * func2(next_x, next_y)) << " | " << std::setprecision(PRECISION / 2) << std::
			setw(10) << func1(next_x, next_y) << " | " << std::setw(10) << func2(next_x, next_y) << std::endl;

		x = next_x;
		y = next_y;
	} while (error > eps);
}

/**
 *  Gradient descent method for solving non linear equation
 */
void gradient_descent_method(double x, double y, std::ostream& ostr) {
	double error;
	int itr = 0;
	int k = 0;
	double vector_of_derivatives[2];
	double lambda = 0.5;
	double alpha;

	ostr << "Newton method" << std::endl << std::endl;
	ostr << "Derivatives matrix:" << std::endl;
	print_derivatives_matrix(ostr);
	ostr << std::endl << std::fixed << "Alpha0 = " << 1 << " Lambda0 = " << 0.5 << std::endl << std::endl;

	//Header
	ostr << "                               Error" << std::
		endl;
	ostr << "  Itr |    x     |     y     | Estimate     |      F1    |      F2    |     FF    | k |   Alpha" <<
		std::endl;

	do {
		itr++;
		// alpha = exp(1); // e = 2.71...
		alpha = 1;
		k = 0;

		double computed_func1 = func1(x, y);
		double computed_func2 = func2(x, y);
		double computed_func_for_minimize = func_for_minimize(x, y);

		vector_of_derivatives[0] = 2 * computed_func1 * func1_derivative_x(x, y) + 2 * computed_func2 *
			func2_derivative_x(x, y);
		vector_of_derivatives[1] = 2 * computed_func1 * func1_derivative_y(x, y) + 2 * computed_func2 *
			func2_derivative_y(x, y);

		while (func_for_minimize(x - alpha * vector_of_derivatives[0], y - alpha * vector_of_derivatives[1]) >=
			computed_func_for_minimize) {
			alpha *= lambda;
			k++;
		}

		double next_x = x - alpha * vector_of_derivatives[0];
		double next_y = y - alpha * vector_of_derivatives[1];

		error = error_estimate_third_norm(x, next_x, y, next_y); 

		ostr << std::fixed << std::setw(4)<< itr << "  | " << std::setprecision(PRECISION) << std::setw(8) << next_x << " | " << std::setw(9) << next_y <<
			" | " << std::setw(10) << std::scientific << sqrt( func1(next_x, next_y) * func1(next_x, next_y) + func2(next_x, next_y) *  func2(next_x, next_y)) << " | " << std::setprecision(PRECISION / 2) << std::
			setw(10) << func1(next_x, next_y) << " | " << std::setw(10) <<  func2(next_x, next_y) << " | " <<
			func_for_minimize(next_x, next_y) << " | " << std::fixed << std::setprecision(PRECISION) << k << " | " << alpha <<
			std::endl;

		x = next_x;
		y = next_y;
	} while (error > eps);
}

/**
 * Print computed jacobian into ostr
 */
void print_computed_jacobian(std::ostream& ostr) {
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
			ostr << std::fixed << std::setw(10) << jacobian[i][j] << " ";
		}

		ostr << std::endl;
	}
}

/**
 * Print equations of the nonlinear system into ostr
 */
void print_nonlinear_system(std::ostream& ostr) {
	switch (variant) {
	case 5:
		ostr << "Fi1(x,y)=sin(x + 0.5) - 1" << std::endl;
		ostr << "Fi2(x,y)=-cos(y - 2)" << std::endl;
		return;
	case 7:
		ostr << "Fi1(x,y)=-sin(x - 1) + 1.3" << std::endl;
		ostr << "Fi2(x,y)=sin(y + 1) + 0.8" << std::endl;
		return;
	case 0:
		ostr << "Fi1(x,y)=-1+sin(y+0.5)" << std::endl;
		ostr << "Fi2(x,y)=-c0s(x-2)" << std::endl;
		return;
	}

	throw std::invalid_argument("Invalid variant");
}

/**
 * Print expressions of the jacobian into ostr
 */
void print_jacobian(std::ostream& ostr) {
	switch (variant) {
	case 5:
		ostr << "cos(x + 0.5) 0" << std::endl;
		ostr << "0            sin(y - 2)" << std::endl;
		return;
	case 7:
		ostr << "-cos(x - 1)  0" << std::endl;
		ostr << "0            cos(y + 1)" << std::endl;
		return;
	case 0:
		ostr << "0.0  ; cos(y+0.5)" << std::endl;
		ostr << "sin(x-2);  0.0 " << std::endl;
		return;
	}

	throw std::invalid_argument("Invalid variant");
}

/**
 * Print expressions of the derivatives matrix into ostr
 */
void print_derivatives_matrix(std::ostream& ostr) {
	switch (variant) {
	case 5:
		ostr << "cos(x + 0.5) -1" << std::endl;
		ostr << "1            sin(y - 2)" << std::endl;
		return;
	case 7:
		ostr << "-cos(x - 1)  1" << std::endl;
		ostr << "1            cos(y + 1)" << std::endl;
		return;
	case 0:
		ostr << " 1.0  ; -cos(y+0.5)" << std::endl;
		ostr << " -sin(x-2);  1.0 " << std::endl;
		return;
	}

	throw std::invalid_argument("Invalid variant");
}
