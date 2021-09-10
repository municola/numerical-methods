#include <iostream>
#include <vector>
#include <cmath>
#include <climits>
#include <float.h>
using namespace std;

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

int main() {
	cout << int_to_bits(20) << endl;
	cout << float_to_bits(4.5) << endl;
}