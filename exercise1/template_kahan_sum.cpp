#include<bits/stdc++.h>

using namespace std;

float naive_sum(vector<float> &v) {
	float sum = 0;
	for(float x : v){
		sum += x;
	}
	return sum;
}

/*
a) If v[0] = 1 (They wrote v[1] = 1 why?), the result is always 1. 
Why? - No clue. 
I thought it would work since float has like 52 digits for the mantissa. 
So 1 + 8.92548 * 10^-8 would be: 
0.10000000000000000.. * 10^1 + 0.892648000000000.. * 10^-8 =
0.100000000892648 * 10^1 (And hence not 1).

If the 1 element is a the end, the summation works perfectly. At least
in my theory. However when i ran it, it wasn't perfectly accurate. Sadly
I coudln't figure out why. 
*/

float sum(vector<float> v) {
	float sum = 0;

	for (int i = 2; i < 1e7; i++) {
		sum += v[i];
	}
	sum += v[0];

	if (v[1] = 1) {
		sum += v[1e7];
		sum += v[1];
	} else {
		sum += v[1];
		sum += v[1e7];
	}

	return sum;
}

/*
c) It's inaccurate again. Both for naive_sum and sum.
naive_sum() results in 1. and sum() in 1.5
*/

float acc_sum(vector<float> &v) {
	float sum = 0;

	//TODO: Point (b)

	return sum;
}

int main() {
	srand(time(0));
	cout << setprecision(15);
	int N = 1e7;
	double corr_sum = 0;
	vector<float> v(N);	
	for (int i = 0; i < N; i++) {
		double x = 1e-8 + (rand() % 10) * 1e-9;
		corr_sum += x;
		v[i] = (float) x;
	}
	corr_sum += 1;
	v[1]=1;

	cout << "naive_sum(v) = " << naive_sum(v) << endl
		 << "      sum(v) = " << sum(v) << endl 
		 << "  acc_sum(v) = " << acc_sum(v) << endl
		 << " correct sum = " << corr_sum << endl;
}
