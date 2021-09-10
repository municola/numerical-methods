#include <iostream>
#include <vector>
#include <cmath>
#include <climits>
#include <float.h>
using namespace std;

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

/*
d)	power and fast_power(double, doulbe) are accurate. fast_power(double, int)
	is not.

e)	Depends on a and b. It gets worse for big b's and/or a's that are very 
	close to a next integer (up rounden) i.e. 3.9999. 
 	It gets better for big a's. i.e 80.
*/

int main() {
	cout << power(2,4) << endl;
	cout << fast_power(2,5) << endl;
	cout << fast_power(2.3, 5) << endl;
}