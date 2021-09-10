#include <iostream>
#include <limits>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Sparse>

using namespace std ;
using namespace Eigen;

int solve_triang(int n){
	MatrixXd A = MatrixXd::Zero(n,n);

	for (int i = 0; i < n; i++) {
		for (int j = i; j < n; j++) {
			A(i,j) = 2;
		}
	}

	FullPivLU<MatrixXd> lu(A);
	cout << lu.matrixLU();

	VectorXd b = 46 *VectorXd::Ones(n);

	VectorXd x = A.fullPivLu().solve(b);
	cout << x << endl;


	return x(n-1);

}

void sample() {
	cout.precision(35);

	double a = 1.0;
	double b = a/9.0;
	double c = b*9.0;

	cout << "b" << b << endl;
	cout << "c" << c << endl;
	if(a==b*9.0) {
		cout << "equal" << endl;
	} else {
		cout << "not equal" << endl;
	}
}

int main() {
	sample();
	return 0;
}
/*
int main() {
	Eigen::Triplet<double> tri(2,3,4);
	MatrixXd A(10,10);
	for (int i = 0; i < 10; i++) {
		A(i,i) = 2;
	}
	for (int i = 5; i < 10; i++) {
		A(i,i) = 3;
	}
	cout << A << endl << endl;
	cout << A.inverse();
	return 0;
}
*/