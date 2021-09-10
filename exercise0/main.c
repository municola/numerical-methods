#include <iostream>
#include <vector>
#include <cmath>
#include <climits>
#include <float.h>
using namespace std;

// Exercise 0.1
void int_to_bits(int x) {
	bool b[32];
	bool a = 1;
	for (int i = 0; i < 32; i++) {
		b[i] = (x & a);
		x = x >> 1;
	}
	
	for (int i = 31; i >= 0; i--) {
		cout << b[i];
	}
	cout << endl; 
}

void float_to_bits(float x) {
	float *a = &x;
	int *b = (int *) a;
	int_to_bits(*b);
}

//---------------------------------------------------------------
// Exercise 0.2
double power (double a, int b) {
	double res = 1;
	for (int i = 0; i < b; i++) {
		res = res * a;
	}
	return res;
}

double fast_power(double a, int b) {
	double res = 1;
	int count = b;
	if (b < 2) {
		return a;
	} else {
		if (fmod(b,2)) {
			res = res * a;
			count--;
		}
		int fp = fast_power(a,count >> 1);
		return res * fp * fp;
	}
}

double fast_power(double a, double b) {
	return exp(b*log(a));
}
// d) power and fast_power(double, doulbe) are accurate. fast_power(double, int)
//	  is not.

// e) Depends on a and b. It gets worse for big b's and/or a's that are very 
//	  close to a next integer (up rounden) i.e. 3.9999. 
// 	  It gets better for big a's. i.e 80.

//----------------------------------------------------------------------
//Exercise 0.3
int factorial(int n) {
	int res = 1;
	for (int i = 1; i <= n; i++) {
		res = res * i;

	}
	return res;
}

double factorial(double n) {
	int a = round(n);
	return factorial(a);
}

// c) for int the result is acccurate till and with 12.
//    for float and double its larger so that DM_MAX >= n!

int main() {
	cout << INT_MAX << endl;
	cout << LONG_MAX << endl;
	cout << DBL_MAX << endl;	
}
