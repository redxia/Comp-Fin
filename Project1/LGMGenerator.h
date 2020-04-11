#ifndef LGMGenerator_
#define LGMGenerator_
#include<cstdlib>
#include<cmath>
#include<vector>
using namespace std;

// LGM constants for the sequence generators (X_n+1 = a * X_n + b) modulo m
const unsigned int a = pow(7, 5); // In LGM Algorithm, a is a constant values: a = 7^5
const unsigned int b = 0; // In LGM Algorithm, b is a constant values: a = 7^5
const unsigned int m = pow(2, 31) - 1; // In LGM Algorithm, b is a constant values: a = 7^5

class LGMGenerator {
public:
	LGMGenerator() { X_0 = 0;}
	// Sets the first number to zero
	void setSeed(const int & x_0);
	int getSeed();


	int LGMNumGenerator(int X);
	vector<double> runif(int n); //Generates the sequence of uniform random numbers
	
	// Generates the sequence of binomial random numbers.
	// p for probability, n for the number of samples, and size is the size of the sequence
	vector<double> rbin(int size, int n, double p); 

	// Generates the exponential distribution sequence.
	vector<double> rexp(int size, double lambda);

	// Generates the normal distribution through Box Muller method
	vector<double> boxMuller(int n);

	// Generates the normal distribution through Polar Marsaglia method
	vector<double> polarMarsaglia(int n);
private:
	int X_0; // Initial Value of the Sequence
	vector<double> unifSeq;
};




#endif 