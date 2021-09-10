	#include <iostream>
#include <algorithm>
#include <vector>
#include <sstream>
#include <string>
#include <complex>

using std::vector;
using std::cout;
using complex = std::complex<double>;

const complex I(0,1); // imaginary unit
const double PI = 3.14159265359;

class Pair {
public:
	int ind;
	complex val;

	Pair(int p, complex v) {
		ind = p;
		val = v;
	}
};

class SparseVec {
public:
	vector<Pair> pairs;
	int len=0;

	SparseVec() {}

	SparseVec(int l) {
		len = l;
	}

	void append(int ind, complex val) {
		// adds val to position ind
		if (ind >= len)
			len = ind+1;

		Pair newPair(ind, val);
		pairs.push_back(newPair);
	}
			
	void cleanup() {
		// sorts w.r.t. indices of pairs and eliminates repeated indices
		std::sort(pairs.begin(), pairs.end(), [] (Pair x, Pair y) { return x.ind < y.ind; });

		vector<Pair> newpairs;
		complex tmp = 0;

		for (int i=0; i<pairs.size(); i++) {
			if (i == pairs.size()-1 || pairs[i+1].ind != pairs[i].ind) {
				tmp += pairs[i].val;
				Pair newpair(pairs[i].ind, tmp);
				newpairs.push_back(newpair);
				tmp = 0;
			}
			else {
				tmp += pairs[i].val;
			}
		}
		pairs = newpairs;
	}

	complex get(int ind) {
		// returns element in position ind
		if (pairs.size() == 0)
			return 0;
		return _get(ind, 0, pairs.size());
	}

	complex _get(int ind, int n1, int n2) {
		if (n2 < n1)
			return 0;
		int m = (n1 + n2)/2;
		if (pairs[m].ind == ind)
			return pairs[m].val;
		if (pairs[m].ind < ind)
			return _get(ind, m+1, n2);
		if (pairs[m].ind > ind)
			return _get(ind, n1, m-1);
	}

	SparseVec conj() const {
		// returns component-wise conjugate
		SparseVec out(len);
		for (auto pair : pairs) {
			out.append(pair.ind, std::conj(pair.val));
		}
		return out;
	}



	static SparseVec cwiseMult(const SparseVec &a, const SparseVec &b) {
		SparseVec out;
		int n = a.len;
		int len = a.pairs.size();
		// TODO: implement component-wise multiplication of SparseVec
		for (int i = 0; i < len; i++) {
			if (a.pairs[i].ind == b.pairs[i].ind){
				Pair newpair(a.paris[i].ind, a.pairs[i].val*b.pairs[i].val);
				out.append(newpair);
			}
		}

		out.cleanup();

		return out;
	}


	static SparseVec conv(const SparseVec &a, const SparseVec &b) {
		SparseVec out;
		int n = a.len;
		int m = b.len;
		out.len = m+n-2;

		n = a.pairs.size();
		m = b.pairs.size();
		// TODO: implement convolution (no conversion to dense structures)
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < m; j++) {
				Pair newpair(a.pairs[i].ind + b.pairs[j].ind,a.pairs[i].val * b.pairs[j].val);
				out.append(newpair);
			}
		}

		out.cleanup();

		return out;
	}

	static SparseVec fft(const SparseVec &x) {
		int n = x.len;
		SparseVec out(n);

		if (n == 1) return x;
		SparseVec even(n/2);
		SparseVec odd(n/2);

		for (int i = 0;i < x.pairs.size(); i++) {
			if (x.pairs[i].ind % 2 == 0){
				even.append(x.pairs[i].ind/2,x.pairs[i].val);
			} else {
				odd.append((x.pairs[i].ind-1)/2,x.pairs[i].val);
			}
		}
		
		SparseVec c1 = fft(even);
		SparseVec c2 = fft(odd);

		c1.cleanup();
		c2.cleanup();

		complex omega = = std::exp(−2∗M_PI/n∗I);

		for (int k = 0; k < n; ++k){
			out.append()
		}


		return out;
	}




	static SparseVec ifft(const SparseVec &x) {
		SparseVec out;

		// TODO: implement inverse fft (no conversion to dense structures)

		return out;
	}

	static SparseVec convfft(SparseVec &a, SparseVec &b) {
		SparseVec out;

		// TODO: implement convolution via fft (no conversion to dense structures)
		
		return out;
	}



	std::string toString() const {
		std::stringstream ss;
		for (auto p : this->pairs) {
			ss << "(" << p.ind << "," << p.val << "),";
		}
		std::string out = ss.str();	
		out+="\n";
		return out;
	}
			

};










/***** TESTING ******/

int main() {

	SparseVec x(5);
	x.append(0,complex(8.2,0));
	x.append(1,complex(1,-2));
	x.append(3,complex(-3,4.66));
	x.append(4,complex(0,4));
	x.cleanup();

	SparseVec y(4);
	y.append(1,complex(5,0));
	y.append(2,complex(1.21,-4));
	y.append(3,complex(4,2.4));
	y.cleanup();

	SparseVec m = SparseVec::cwiseMult(x,y);
	m.cleanup();
	cout << "TESTS. Correct componentwise multiplication between x and y: (1,(5,-10)),(3,(-23.184,11.44)),\n";
	cout << "cwiseMult(x,y) = " << m.toString();

	SparseVec c = SparseVec::conv(x,y);
	c.cleanup();
	cout << "Correct exact discrete convolution between x and y: (1,(41,0)),(2,(14.922,-42.8)),(3,(26.01,13.26)),(4,(-6.2,17.7)),(5,(15.01,37.6386)),(6,(-7.184,16.28)),(7,(-9.6,16)),\n";
	cout << "conv(x,y) = " << c.toString();
	SparseVec cf = SparseVec::convfft(x,y);
	cf.cleanup();
	cout << "convfft(x,y) = " << cf.toString();
}

