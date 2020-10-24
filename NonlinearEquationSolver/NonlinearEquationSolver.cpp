#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include "Functions.h"
using namespace std;

const bool PRINT_TO_FILE = true;

int main() {
	ofstream fout("../" + to_string(variant) + "_output" + ".txt");
	ostream& out = PRINT_TO_FILE ? fout : cout;

	out << "Variant " << variant << endl;
	out << "b" << endl;

	double x;
	double y;

	set_x0_y0(x, y);
	out << "x0 = " << x << ", y0 = " << y << endl;

	simple_iteration_method(x, y, out);
	out << endl;

	newton_method(x, y, out);
	out << endl;

	gradient_descent_method(x, y, out);
	out << endl;

	fout.close();
    return 0;
}