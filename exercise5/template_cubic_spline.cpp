// to compile run: g++ -std=gnu++11 cubic_spline.cpp -lmgl
#include <iostream>
#include <eigen3/Eigen/Dense>
#include <mgl2/mgl.h>

using namespace Eigen;

MatrixXd cubicSpline(const VectorXd &T, const VectorXd &Y) {
	// returns the matrix representing the spline interpolating the data
	// with abscissae T and ordinatae Y. Each column represents the coefficients
	// of the cubic polynomial on a subinterval.
	// Assumes T is sorted, has no repeated elements and T.size() == Y.size().

	int n = T.size() - 1; // T and Y have length n+1

	// TODO: build the spline matrix with polynomials' coefficients
	MatrixXd spline(4, n);
	n++; // SORRY. messed up. thought n = length n+1;

	// Vector h
	VectorXd h;
	for (int i = 1; i < n; i++) {
		h(i) = T(i)-T(i-1); // h(1) = t1 - t0, h(0) = buffer
	}

	// Matrix for LSE to determine the rhos
	MatrixXd A;
	// Diagonal
	for (int i = 0; i < n-1; i++) {
		A(i,i) = (h(i+1)+h(i+2))/3;
	}
	// left Diagonal and right Diagonal
	for (int i = 0; i < n-2; i++) {
		A(i+1,i) = h(i+2)/6;
		A(i,i+1) = h(i+2)/6;
	}

	// Vector r for LSE to determine the rhos (length n-1)
	VectorXd r;
	for (int i = 1; i < n-1; i++) {
		r(i-1) = (Y(i+1)-Y(i))/h(i+1)-(Y(i)-Y(i-1))/h(i); // formular but adpated by -1 bc we start at 0;
	}

	// Solve the LSE. Ax=r, rho = x
	VectorXd x;
	x = A.lu().solve(r);

	// now we have the rhos wo we can construct eh Matrix spline
	// do the a's
	for (int i = 0; i < n-1; i++) {
		spline(0,i) = Y(i); 
	}
	// do the b's
	for (int i = 0; i < n-1; i++) {
		spline(1,i) = (Y(i+1)-Y(i))/h(i+1) - (h(i+1)*(2*x(i)+x(i+1)))/6;
	}
	// do the c's
	for (int i = 0; i < n-1; i++) {
		spline(2,i) = x(i)/2;
	}
	// do the d's
	for (int i = 0; i < n-1; i++) {
		spline(3,i) = (x(i+1)-x(i))/(6*h(i+1));
	}

	return spline;
}

VectorXd evalCubicSpline(const MatrixXd &S, const VectorXd &T, const VectorXd &evalT) {
	// Returns the values of the spline S calculated in the points X.
	// Assumes T is sorted, with no repetetions.

	int n = evalT.size();
	VectorXd out(n);

	// TODO: fill out
	// Find the right intervall
	int last = 0;
	for (int i = 0; i<n;i++) { // For every evalT
		for (int j = 0; j < T.size()-1;j++) {
			if (evalT(i)>T(i) && evalT(i)<T(i+1)) {
				double x = evalT(i)-T(j);
				out(i) = S(0,j) + x*(S(1,j)+x*(S(2,j) + x*(S(3,j))));
				break;
			}
		}
	}


	return out;
}

int main() {
	// tests
	VectorXd T(9);
	VectorXd Y(9);
	T << 0, 0.4802, 0.7634, 1, 1.232, 1.407, 1.585, 1.879, 2;
	Y << 0., 0.338, 0.7456, 0, -1.234, 0 , 1.62, -2.123, 0;

	int len = 1 << 9;
	VectorXd evalT = VectorXd::LinSpaced(len, T(0), T(T.size()-1));

	VectorXd evalSpline = evalCubicSpline(cubicSpline(T, Y), T, evalT);

 	mglData datx, daty;
	datx.Link(evalT.data(), len);
	daty.Link(evalSpline.data(), len);
	mglGraph gr;
	gr.SetRanges(0, 2, -3, 3);
	gr.Plot(datx, daty, "0");
	gr.WriteFrame("spline.eps");
}
	
