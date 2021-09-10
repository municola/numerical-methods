#include <iostream>
#include <algorithm>
#include <vector>
#include <sstream>
#include <string>
#include <complex>

using namespace std;

const complex<double> I(0,1); // imaginary unit
const double PI = 3.14159265359;

template<class T> struct duplet {
	int ind;
	T val;

	duplet(int p, T v) {
		ind = p;
		val = v;
	}
};

template<class T> struct sparse_vec {
	double tol = 1e-6;
	vector<duplet<T> > duplets;
	int len=0;

	sparse_vec(int l) {
		len = l;
	}

	void append(int ind, T val) {
		if (std::abs(val) >= tol) {
			newDublet = duplet(ind,val);
			duplets.push_back(newDublet);
		} 
	}
			
	void cleanup() {
		// first sort
		std::sort()::sort(duplets.begin(), duplets.end(), 
			[](duplet<T> a, duplet<T> b) {
				return a.ind < b.ind;
			}
		);

		// destroy dulicats AND
                vector<duplet<T> > newduplets;
                T tmp = 0;
                for (int i=0; i<duplets.size(); i++) {
                        if (i == duplets.size()-1 || duplets[i+1].ind != duplets[i].ind) {
                                tmp += duplets[i].val;
                                if (abs(tmp) >= tol) {
                                        duplet<T> newduplet(duplets[i].ind, tmp);
                                        newduplets.push_back(newduplet);
                                }
                                tmp = 0;
                        }
                        else {
                                tmp += duplets[i].val;
                        }
                }
                duplets = newduplets;
	}

	T get_val(int ind) const {
		//TODO
		return 0;
	}

	static sparse_vec cwise_mult(const sparse_vec &a, const sparse_vec &b) {
		sparse_vec out(max(a.len,b.len));
		//TODO
		return out;
	}

	static sparse_vec conv(const sparse_vec &a, const sparse_vec &b) {
		sparse_vec out(a.len + b.len - 1);
		//TODO
		return out;
	}

	static sparse_vec fft(const sparse_vec &x) {
		int n = x.len;
		sparse_vec tot(n);
		//TODO
		return tot;
	}

	static sparse_vec ifft(const sparse_vec &x) {
		double n = x.len;
		sparse_vec out(n);
		//TODO
		return out;
	}

	static sparse_vec conv_fft(sparse_vec a, sparse_vec b) {
		//TODO
		return a;
	}

	std::string to_string() const {
		std::stringstream ss;
		for (auto p : this->duplets) {
			ss << "(" << p.ind << "," << p.val << "),";
		}
		ss << "\n";
		std::string out = ss.str();	
		return out;
	}
			

};




/***** TESTING ******/

int main() {

	sparse_vec<complex<double> > x(5);
	x.append(0,complex<double>(8.2,0));
	x.append(1,complex<double>(1,-2));
	x.append(3,complex<double>(-3,4.66));
	x.append(4,complex<double>(0,4));
	x.cleanup();

	sparse_vec<complex<double> > y(4);
	y.append(1,complex<double>(5,0));
	y.append(2,complex<double>(1.21,-4));
	y.append(3,complex<double>(4,2.4));
	y.cleanup();

	auto m = sparse_vec<complex<double> >::cwise_mult(x,y);
	m.cleanup();
	cout << "TESTS. Correct componentwise multiplication between x and y: (1,(5,-10)),(3,(-23.184,11.44)),\n";
	cout << "cwise_mult(x,y) = " << m.to_string();

	auto c = sparse_vec<complex<double> >::conv(x,y);
	c.cleanup();
	cout << "Correct exact discrete convolution between x and y: (1,(41,0)),(2,(14.922,-42.8)),(3,(26.01,13.26)),(4,(-6.2,17.7)),(5,(15.01,37.6386)),(6,(-7.184,16.28)),(7,(-9.6,16)),\n";
	cout << "conv(x,y) = " << c.to_string();
	auto cf = sparse_vec<complex<double> >::conv_fft(x,y);
	cf.cleanup();
	cout << "conv_fft(x,y) = " << cf.to_string();
}


