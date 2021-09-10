#include <iostream>

#include <eigen3/Eigen/Dense>

using namespace Eigen;
using namespace std;

struct Newton {
	Newton(const Eigen::VectorXd &x) : _x(x), _a(x.size()) { }
	void Interpolate(const Eigen::VectorXd &y);
	double operator()(double x) const;

private:
	Eigen::VectorXd _x;	// nodes
	Eigen::VectorXd _a;	// coefficients
};

// Compute the coefficients in the Newton basis.
void Newton::Interpolate(const Eigen::VectorXd &y) {
	// TODO: Task (a)
	int n = y.size();

	// Construct Matrix A
	MatrixXd A(n,n);
	VectorXd c = VectorXd::Ones(n);
	A.col(0) = c;
	for (int i = 1;i<n;i++) { // i = cols, j = rows
		for (int j = 0; j <n; j++) {
			A(j,i) = A(j,i-1)*(_x[j]-_x[i-1]);
		}
	} 

	VectorXd a(n);

	_a = A.triangularView<Lower>().solve(y);

}

// Evaluate the interpolant at x.
double Newton::operator()(double x) const {
	int n = _x.size();
	double result = _a[n-1];
	for (int i = n -2; i >= 0; i--) {
		result *= (x-_x[i]);
		result += _a[i];
	}


	return result; // dummy
}

struct Lagrange {
	Lagrange(const Eigen::VectorXd &x);
	void Interpolate(const Eigen::VectorXd &y) { _y = y; }
	double operator()(double x) const;

private:
	Eigen::VectorXd _x;	// nodes
	Eigen::VectorXd _l;	// weights
	Eigen::VectorXd _y;	// coefficients
};

// Compute the weights l for given nodes x.
Lagrange::Lagrange(const Eigen::VectorXd &x) : _x(x), _l(x.size()), _y(x.size()) {
	// TODO: Task (c)
	int n = x.size();

	for (int i = 0; i < n; i++) {
		_l[i] = 1;
		for (int j = 0; j < n; j++) {
			if (i!=j) {
				_l[i] *= 1/(_x[i]-_x[j]);
			}
		}
	}
	// compute _l
	// ...
}

// Evaluate the interpolant at x.
double Lagrange::operator()(double x) const {
	// TODO: Task (d)
	int n = _x.size();
	// computes the w(x)
	double w = 1;
	for (int i = 0; i < n; i++) {
		w *= (x-_x[i]);
	}

	// Calculate the polynomial
	double result = 0;
	for (int i = 0; i<n; i++) {
		result += _y[i]*w*_l[i]/(x-_x[i]);
	}
	return result; // dummy
}

// Runge function
Eigen::VectorXd r(const Eigen::VectorXd &x) {
	return (1.0 / (1.0 + 25.0 * x.array() * x.array())).matrix();
}

int main() {
	int n = 5;
	Eigen::VectorXd x;
	x.setLinSpaced(5, -1.0, 1.0);
	Eigen::VectorXd y = r(x);

	Newton p(x);
	p.Interpolate(y); // correct result: p._a = [0.0384615, 0.198939, 1.5252, -3.31565, 3.31565]

	Lagrange q(x);    // correct result: p._l = [0.666667, -2.66667, 4, -2.66667, 0.666667]
	q.Interpolate(y);

	// Compute difference of p and q.
	int m = 22;
	double offset = 0.08333333333;
	x.setLinSpaced(m, -1.0 + offset, 1.0 - offset);
	double norm2 = .0;
	for (int i = 0; i < m; ++i) {
		double d = p(x(i)) - q(x(i));
		norm2 += d * d;
	}

	// By uniquenss of the interpolation polynomial, we expect p = q.
	std::cout << "This number should be close to zero: " << norm2 << std::endl;

	return 0;
}

