#include <iostream>
#include <vector>
#include <cmath>
#include <climits>
#include <float.h>
using namespace std;

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

/*
c) for int the result is acccurate till and with 12.
   for float and double its larger.
*/

int main() {
	cout << factorial(5) << endl;
	cout << factorial(5.3) << endl;
}
