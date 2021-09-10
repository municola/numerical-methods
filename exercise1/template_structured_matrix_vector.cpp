// compile with: g++ template_structured_matrix_vector.cpp -I/usr/include/eigen3 -lmgl
// -std=c++11

#include <chrono>
#include <iostream>
#include <iomanip>
//#include <limits>
//#include <ratio>
#include <vector>

#include <Eigen/Dense>
#include <mgl2/mgl.h>

using namespace Eigen;
using namespace std;

typedef std::chrono::time_point<std::chrono::high_resolution_clock> Timepoint;
typedef std::chrono::high_resolution_clock Time;
typedef std::chrono::duration<double, std::ratio<1>> Duration;

/* \brief compute $\mathbf{A}\mathbf{x}$

 * \mathbf{A} is defined by $(\mathbf{A})_{i,j} := \min {i,j}$
 * \param[in] x vector x for computation of A*x = y
 * \param[out] y = A*x
 */
void multAminSlow(const VectorXd & x, VectorXd & y) {
    unsigned int n = x.size();

    VectorXd one = VectorXd::Ones(n);
    VectorXd linsp = VectorXd::LinSpaced(n,1,n);
    y = ( ( one * linsp.transpose() )
          .cwiseMin( linsp * one.transpose()) ) * x;
}

/* \brief compute $\mathbf{A}\mathbf{x}$
 * \mathbf{A} is defined by $(\mathbf{A})_{i,j} := \min {i,j}$
 * Instead of a "Matlab style" construcion of the product,
 * we use simple loops.
 * \param[in] x vector x for computation of A*x = y
 * \param[out] y = A*x
 */
void multAminLoops(const VectorXd & x, VectorXd & y) {
    unsigned int n = x.size();
    MatrixXd A(n,n);

    for(unsigned int i = 0; i < n; ++i) {
        for(unsigned int j = 0; j < n; ++j) {
            A(i,j) = std::min(i+1,j+1);
        }
    }
    y = A * x;
}

/* \brief compute $\mathbf{A}\mathbf{x}$
 * This function has optimal complexity.
 * \mathbf{A} is defined by $(\mathbf{A})_{i,j} := \min {i,j}$
 * \param[in] x vector x for computation of A*x = y
 * \param[out] y = A*x
 */
void multAmin(const VectorXd & x, VectorXd & y) {
    unsigned int n = x.size();
    y = VectorXd::Zero(n);
    VectorXd z = VectorXd::Zero(n);
    z(n-1) = x(n-1);
    for (int i = n-2; i >= 0; i--) {
      z(i) = z(i+1) + x(i);
    }
    y(0) = z(0);
    for (int i = 1; i < n; i++) {
      y(i) = z(i) + y(i-1);
    }
}

int main(void) {
    // Testing correctness of the code
    unsigned int M = 10;
    VectorXd xa = VectorXd::Random(M);
    VectorXd ys, yf;

    multAmin(xa, yf);
    multAminSlow(xa, ys);
    // Error should be small
    std::cout << "||ys-yf|| = " << (ys - yf).norm() << std::endl;


    unsigned int nLevels = 9;
	  unsigned int *n = new unsigned int[nLevels];
	  double *minTime = new double[nLevels];
	  double *minTimeLoops = new double[nLevels];
	  double *minTimeEff = new double[nLevels];

	  n[0] = 4;
	  for (unsigned int i=1; i<nLevels; i++)
		  n[i] = 2*n[i-1];

	   //TODO: Point (c)
    unsigned int N = 16;
    for (int i = 0; i < 7; i++) {
      double minTime = numeric_limits<double>::infinity();
      for (int j = 0; j < 10; j++) {
        VectorXd xa = VectorXd::Random(N);
        VectorXd ys, yf;
        Timepoint start = Time::now();

        multAmin(xa,ys);

        Timepoint end = Time::now();
        Duration d = end - start;
        minTime = min(minTime,d.count());
      }
      cout << i << ": " << minTime << std::scientific << std::setprecision(3) << endl;
      N = N*2;
    }
    
    // Plotting with MathGL
    double nMgl[nLevels];
    double ref1[nLevels], ref2[nLevels];
    for (int i=0; i<nLevels; i++) {
    	nMgl[i] = n[i];
    	ref1[i] = 1e-8*pow(n[i],2);
    	ref2[i] = 1e-7*n[i];
    }
    
    mglData matSize;
    matSize.Link(nMgl, nLevels);
    
    mglData data1, data2;
    mglData dataRef1, dataRef2;
  	data1.Link(minTime, nLevels);
  	data2.Link(minTimeEff, nLevels);
  	dataRef1.Link(ref1,nLevels);
  	dataRef2.Link(ref2,nLevels);
  	
  	mglGraph *gr = new mglGraph;
    gr->Title("Runtime of multAmin");
  	gr->SetRanges(n[0],n[0]*pow(2,nLevels-1),1e-6,1e+1);  gr->SetFunc("lg(x)","lg(y)");
  	gr->Axis();
  	gr->Plot(matSize,data1,"k +"); gr->AddLegend("slow","k +");
  	gr->Plot(matSize,data2,"r +"); gr->AddLegend("efficient","r +");
  	gr->Plot(matSize,dataRef1,"k"); gr->AddLegend("O(n^2)","k");
  	gr->Plot(matSize,dataRef2,"r"); gr->AddLegend("O(n)","r");
  	gr->Label('x',"Matrix size [n]",0);
  	gr->Label('y', "Runtime [s]",0);
    gr->Legend(2);
	  gr->WriteFrame("multAmin_comparison.eps");


    // The following code is just for demonstration purposes.
    // Build Matrix B with dimension 10x10
    unsigned int nn = 10;
    MatrixXd B = MatrixXd::Zero(nn,nn);
    for(unsigned int i = 0; i < nn; ++i) {
        B(i,i) = 2;
        if(i < nn-1) B(i+1,i) = -1;
        if(i > 0) B(i-1,i) = -1;
    }
    B(nn-1,nn-1) = 1;

    //TODO: Point (e)
    MatrixXd A(10,10);
    for(unsigned int i = 0; i < 10; ++i) {
        for(unsigned int j = 0; j < 10; ++j) {
            A(i,j) = std::min(i+1,j+1);
        }
    }

    for (int i = 0; i < 10; i++) {
      VectorXd e = VectorXd::Zero(10);
      e(i) = 1;
      cout << A*B*e << endl;
      cout << " " << endl;
    }

    cout << A*B;

}

/*
a) O(n²)
e) A*B = I
*/