#include <chrono>
#include <functional>
#include <iostream>
#include <vector>

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseLU>

using namespace std;

std::vector<Eigen::Triplet<double>> MakeTripletList(int n) {
	typedef Eigen::Triplet<double> T;

	int nnz = 3 * n - 2;
	std::vector<Eigen::Triplet<double>> tripletList(nnz);

	for (int i = 0; i < n-1; i++) {
		tripletList.push_back(T(i,i,2.0));
		tripletList.push_back(T(i+1,i,-1.0));
		tripletList.push_back(T(i,i+1,-1.0));
	}

	tripletList.push_back(T(n-1,n-1,2.0));

	return tripletList;
}

double Runtime(const std::function<void(void)> &f) {

	std::chrono::time_point<std::chrono::high_resolution_clock> start;
	std::chrono::time_point<std::chrono::high_resolution_clock> end;

	// warmup
	f();

	// Execute 10 times and take average
	double duration = 0.0;
	int counter = 10;

	for (int i = 0; i < counter; i++) {
		 start = std::chrono::high_resolution_clock::now();
		 f();
		 end = std::chrono::high_resolution_clock::now();
		 duration = duration + (double)(end - start).count();
	}

    return duration/counter;	// dummy return value
}

template <class T>
std::ostream & operator<< (std::ostream &os, const std::vector<T> &v) {
	os << "[";
	if (!v.empty()) {
		os << v[0];
		for (int i = 1; i < v.size(); ++i) os << ", " << v[i];
	}
    os << "]";

    return os;
}

int main() {
	// print small example of the tridiagonal matrix
	int m = 4;
	std::vector<Eigen::Triplet<double>> tripletList = MakeTripletList(m);
	Eigen::SparseMatrix<double> S_(m, m);
	S_.setFromTriplets(tripletList.begin(), tripletList.end());
	std::cout << "If n = " << m << ", then T equals" << std::endl;
	std::cout << Eigen::MatrixXd(S_) << std::endl;

	// matrix sizes for benchmark
	std::vector<int> N = {64, 128, 256};
	std::cout << "LU decomposition of T, where n = " << N << std::endl;

	// set up variables for runtime measurement
	std::vector<double> runtimeSparse;
	std::vector<double> runtimeDense;

	for (int n : N) {
		tripletList = MakeTripletList(n);

		// sparse LU decomposition
	    Eigen::SparseMatrix<double> S(n, n);
	    S.setFromTriplets(tripletList.begin(), tripletList.end());
	    Eigen::SparseLU<Eigen::SparseMatrix<double>> sparseLU;
	    std::function<void(void)> SparseLU = [&S, &sparseLU] () {
	            sparseLU.compute(S);
	    };

	    // dense LU decomposition
	    Eigen::MatrixXd D(S);
	    Eigen::FullPivLU<Eigen::MatrixXd> denseLU(n, n);
	    std::function<void(void)> DenseLU = [&D, &denseLU] () {
	            denseLU.compute(D);
	    };

	    // benchmark
	    double timeSparse = Runtime(SparseLU);
	    runtimeSparse.push_back(timeSparse);
	    double timeDense = Runtime(DenseLU);
	    runtimeDense.push_back(timeDense);
	}

	std::cout << "Runtime in seconds using storage format..." << std::endl;
	std::cout << "...sparse: " << runtimeSparse << std::endl;
	std::cout << "...dense:  " << runtimeDense << std::endl;

	return 0;
}
