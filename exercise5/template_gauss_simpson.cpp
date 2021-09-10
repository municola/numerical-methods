#include <cmath>
#include <iostream>
#include <eigen3/Eigen/Dense>
#include <functional>
#include <vector>

using namespace Eigen;
using namespace std;

double GaussLegendre5(const std::function<double(double)> &f, double a, double b) {
	VectorXd c(5);
	VectorXd w(5);

	c <<  -0.90617984593, -0.53846931010, 0.0, 0.53846931010, 0.90617984593;
	w << 0.23692688505, 0.47862867049, 0.56888888888, 0.47862867049, 0.23692688505;

	// Adapt c and w to the right intervall
	for (int i = 0; i<c.size(); i++) {
		c(i) = (0.5*(1-c(i))*a + 0.5*(1+c(i))*b);
		w(i) = 0.5*(b-a)*w(i);
	}

	// Evaluate
	double result = 0;
	for (int i = 0; i < c.size(); i++) {
		result += w(i)*f(c(i));
	}

	return result; // dummy
}

double CompositeSimpson(const std::function<double(double)> &f, const std::vector<double> &x) {

	// Compute f(x) and store it in vecote fx --fx(i) stores f(xj)
	VectorXd fx(x.size());
	for (int i = 0; i < x.size(); i++) {
		fx(i) = f(x(i));
	}

	// Compute f(0.5*(xj-xj-1)) and store it in vector fxx --fxx(1) stores f(x1-x0)
	VectorXd fxx(x.size());
	for (int i = 1; i < x.size(); i++) {
		fxx(i) = f(0.5*(x(i)-x(i-1)))
	}

	// Compute f(a) and f(b)
	VectorXd fab(2);
	fab(0) = f(a);
	fab(1) = f(b);

	result = 0;

	// First summation before the sum sign
	result += 1.0/6.0*(x(1)-x(0))*fab(0);

	// Summation for the first summ sign
	for (int i = 1; i < x.size()-1; i++) {
		result += 1.0/6.0*(x(i+1)-x(i-1))*fx(i);
	}

	// Summation of the second summ sign
	for (int i = 1; i < x.size(); i++) {
		result += 2.0/3.0*(x(i)-x(i-1))*fxx(i)
	}

	// Last summation. that one after the last sum sign
	result += 1.0/6.0*(x(x.size()-1)-x.size()-2)*fab(1);
	

	return result;
}

std::vector<double> LinSpace(int n, double a, double b) {
	std::vector<double> x(n);
	double d = (b - a) / (n - 1);
	for (int i = 0; i < n; ++i) {
		x[i] = i * d + a;
	}
	return x;
}

double f(double x) {
	return std::sqrt(x);
}

double F(double x) {
	double y = std::sqrt(x);
	return 2.0 / 3.0 * y * y * y;
}

int main() {
	int n = 5;
	double a = 0.0;
	double b = 1.0;

	std::cout << "Gauss-Legendre: " << GaussLegendre5(f, a, b) << std::endl;
	std::cout << "Simpson: " << CompositeSimpson(f, LinSpace(n, a, b)) << std::endl;
	std::cout << "Exact value: " << F(b) - F(a) << std::endl;

	return 0;
}
