#include<bits/stdc++.h>
#include<eigen3/Eigen/Dense>

using namespace std;
using namespace Eigen;

MatrixXf mult(const MatrixXf &A, const MatrixXf &B) {
	int N = A.rows();
	MatrixXf C = MatrixXf::Zero(N, N);

	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			for (int  k = 0; k < N; k++) {
				C(i,j) += A(i,k)*B(k,j); 
			}
		}
	}

	return C;
}

MatrixXf mult_rec(const MatrixXf &A, const MatrixXf &B) {
	int N = A.rows();
	MatrixXf C(N, N);
	if (N == 2) {
		C(0,0) = A(0,0)*B(0,0) + A(0,1)*B(1,0);
		C(0,1) = A(0,0)*B(0,1) + A(0,1)*B(1,1);
		C(1,0) = A(1,0)*B(0,0) + A(1,1)*B(1,0);
		C(1,1) = A(1,0)*B(0,1) + A(1,1)*B(1,1);
	} else {
		C.block(0,0,N/2,N/2) = mult_rec(A.block(0,0,N/2,N/2),B.block(0,0,N/2,N/2)) + 
			mult_rec(A.block(0,N/2,N/2,N/2),B.block(N/2,0,N/2,N/2));
		C.block(0,N/2,N/2,N/2) = mult_rec(A.block(0,0,N/2,N/2),B.block(0,N/2,N/2,N/2)) + 
			mult_rec(A.block(0,N/2,N/2,N/2),B.block(N/2,N/2,N/2,N/2));
		C.block(N/2,0,N/2,N/2) = mult_rec(A.block(N/2,0,N/2,N/2),B.block(0,0,N/2,N/2)) + 
			mult_rec(A.block(N/2,N/2,N/2,N/2),B.block(N/2,0,N/2,N/2));
		C.block(N/2,N/2,N/2,N/2) = mult_rec(A.block(N/2,0,N/2,N/2),B.block(0,N/2,N/2,N/2)) + 
			mult_rec(A.block(N/2,N/2,N/2,N/2),B.block(N/2,N/2,N/2,N/2));
	}
	return C;
}
	

MatrixXf strassen(const MatrixXf &A, const MatrixXf &B) {
	int N = A.rows();
	MatrixXf C(N, N);

	//TODO: Point (e)
	if (N == 1) {
		C = A*B;
	} else {
		#define b11 B.block(0,0,N/2,N/2)
		#define b12 B.block(0,N/2,N/2,N/2)
		#define b21 B.block(N/2,0,N/2,N/2)
		#define b22 B.block(N/2,N/2,N/2,N/2)

		#define a11 A.block(0,0,N/2,N/2)
		#define a12 A.block(0,N/2,N/2,N/2)
		#define a21 A.block(N/2,0,N/2,N/2)
		#define a22 A.block(N/2,N/2,N/2,N/2)

		MatrixXf M1 = strassen((a11+a22),(b11+b22));
		MatrixXf M2 = strassen((a21+a22),b11);
		MatrixXf M3 = strassen(a11,(b12-b22));
		MatrixXf M4 = strassen(a22,(b21-b11));
		MatrixXf M5 = strassen((a11+a12),b22);
		MatrixXf M6 = strassen((a21-a11),(b11+b12));
		MatrixXf M7 = strassen((a12-a22),(b21+b22));

		C.block(0,0,N/2,N/2) = M1 + M4 - M5 + M7;
		C.block(0,N/2,N/2,N/2) = M3 + M5;
		C.block(N/2,0,N/2,N/2) = M2 + M4;
		C.block(N/2,N/2,N/2,N/2) = M1 - M2 + M3 + M6;
	}

	return C;
}

int main() {
	srand(time(0));
	cout << setprecision(6) << setfill(' ');

	for (int i = 1; i < 9; i++) {
		int N = 1 << i;		
		cout << "Matrix size = " << N << endl;
		MatrixXd AA = MatrixXd::Random(N, N);
		MatrixXd BB = MatrixXd::Random(N, N);
		MatrixXd ans = AA*BB;
		MatrixXf A = AA.cast<float>();
		MatrixXf B = BB.cast<float>();

		auto start = std::chrono::steady_clock::now();
		MatrixXf W = mult(A, B);
		auto finish = std::chrono::steady_clock::now();
		cout << setw(24) << " " <<  setw(15) 
			<< "Time (s)" << setw(20) << "Error (l2-norm)"  << endl;
		cout << setw(24) << "Naive iterative "<< setw(15) 
			<< std::chrono::duration_cast<std::chrono::duration<double> > 
			(finish - start).count() << setw(20) << (W.cast<double>() - ans).norm() << endl;
	
		start = std::chrono::steady_clock::now();
		MatrixXf X = mult_rec(A, B);
		finish = std::chrono::steady_clock::now();
		cout << setw(24) << "Naive recursive "<< setw(15) 
			<< std::chrono::duration_cast<std::chrono::duration<double> > 
			(finish - start).count() << setw(20) << (X.cast<double>() - ans).norm() << endl;
		start = std::chrono::steady_clock::now();
		MatrixXf Y = strassen(A, B);
		finish = std::chrono::steady_clock::now();
		cout << setw(24) << "Strassen recursive " << setw(15) 
			<< std::chrono::duration_cast<std::chrono::duration<double> >
			(finish - start).count() << setw(20) << (Y.cast<double>() - ans).norm() << endl;
		start = std::chrono::steady_clock::now();	
		MatrixXf Z = A*B;
		finish = std::chrono::steady_clock::now();
		cout << setw(24) << "Eigen built-in "<< setw(15) 
			<< std::chrono::duration_cast<std::chrono::duration<double> >
			(finish - start).count() << setw(20) << (Z.cast<double>() - ans).norm() << "\n\n\n";
	}
}

/*
b) O(n³)
*/