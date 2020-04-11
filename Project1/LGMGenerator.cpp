/*
 * Author: Redmond Xia
 * Date: 04/7/2020
 * Details: Here we implement the Lewis, Goodman, and Miller algorithms for Random Numbers
 */

#include "LGMGenerator.h"
#include<iostream>
#include<cmath>

const double PI = 3.159265358979;

void LGMGenerator::setSeed(const int & x_0) {
	X_0 = x_0;
}

int LGMGenerator::getSeed() {
	return X_0;
}

int LGMGenerator::LGMNumGenerator(int X) {
	return (a * X + b) % m; //a,b,m are global constant found in LGMGenerator.h
}

vector<double> LGMGenerator::runif(int n) {
	vector<double> X(n);
	X[0] = getSeed();
	// Generates the X sequences first
	for (int i = 1; i < X.size(); i++) {
		X[i] = LGMNumGenerator(X[i-1]);
	}

	// Adjust the domain to [0,1]
	for (int i = 0; i < X.size(); i++) {
		X[i] /= double(m);
	}
	return X;
}

vector<double> LGMGenerator::rbin(int size, int n, double p) {
	vector<double> X(n);
	vector<double> binomial(size);
	
	for (int i = 0; i < size; i++) {
		int sum = 0;
		setSeed(i+1);
		X = runif(n);
		// Summing the Bernoulli sequence
		for (int j = 0; j < n; j++) {
			if (X[j] <= p) {
				sum += 1;
			}
		}
		binomial[i] = sum;
	}
	return binomial;
}

vector<double> LGMGenerator::rexp(int n, double lambda) {
	vector<double> X(n);
	X = runif(n);
	vector<double> arr(n);
	for (int i = 0; i < X.size(); i++) {
		arr[i] = -lambda * log(X[i]);
	}
	return arr;
}

vector<double> LGMGenerator::boxMuller(int n) {
	vector<double> X(n);
	X = runif(n);

	// Stores the normal distribution, in sequence as Z1,Z2,Z1,Z2 since can only return one vector
	vector<double> arr(n);
	for (int i = 0; i < arr.size(); i += 2) {
		// These formula come from the Box-Muller Method
		arr[i] = sqrt(-2 * log(X[i])) * cos(2 * PI * X[i + 1]);
		arr[i+1] = sqrt(-2 * log(X[i])) * sin(2 * PI * X[i + 1]);
	}
	return arr;
}

vector<double> LGMGenerator::polarMarsaglia(int n) {
	vector<double> X(n);
	X = runif(n);

	// Stores the normal distribution, in sequence as Z1,Z2,Z1,Z2 since can only return one vector
	vector<double> arr(n);
	int counter = 0;
	double v1 = 0;
	double v2 = 0;
	double w = 0;
	for (int i = 0; i < arr.size(); i += 2) {
		v1 = 2 * X[i] - 1;
		v2 = 2 * X[i + 1] - 1;
		w = v1 * v1 + v2 * v2;

		if (w <= 1.0) {
			arr[counter] = v1 * sqrt(-2 * log(w) / w);
			arr[counter+1] = v2 * sqrt(-2 * log(w) / w);
			counter += 2;
		}
	}
	return arr;
}