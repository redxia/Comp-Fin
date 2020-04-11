/*
 * Author: Redmond Xia
 * Date: 04/07/2020
 *
 *
 *
 *
 */

#include<iostream>
#include<string>
#include<cmath>
#include<cstdlib>
#include "LGMGenerator.h"
#include<vector>
#include<random>
#include<fstream>
#include<chrono>

using namespace std;
using namespace std::chrono;

// Helper Functions
void printVector(vector<double> X);
template <class T>
double mean(vector<T> X);
template <class T>
double stdDev(vector<T> X);

//Creating a txt to plot in R
template <class T>
void writeToTXT(vector<T> arr, string filename);

// Returns the sequence of independent bernouilli trials based on the uniform values
vector<int> indepBern(vector<double> X, int n);
int probBin(vector<double> X);


int main() {
	class LGMGenerator Generator;
	Generator.setSeed(5); // Set the initial value of the sequence

	cout << "============================== Question 1  ==============================" << endl;
	cout << "============================== Question 1a ==============================" << endl;
	int n = 10000; // Size of the sequence
	vector<double> X_n(n);
	X_n = Generator.runif(n); // Generates the uniform distribution sequence
	cout << "The empirical mean with X_0 = " << Generator.getSeed() << ": " << mean(X_n) << endl;
	cout << "The empirical standard deviation with X_0 = " << Generator.getSeed() << ": " << stdDev(X_n) << endl;
	cout << "=========================================================================" << endl << endl;

	cout << "============================== Question 1b ==============================" << endl;
	srand(0);
	vector<double> unif(n);
	default_random_engine libGenerator;
	uniform_real_distribution<double> distribution(0.0, 1.0);
	for (int i = 0; i < unif.size(); i++) {
		unif[i] = distribution(libGenerator);
	}
	cout << "The built in library uniform distribution mean: " << mean(unif) << endl;
	cout << "The built in library uniform distribution standard deviation: " << stdDev(unif) << endl;
	cout << "=========================================================================" << endl << endl;

	cout << "============================== Question 1c ==============================" << endl;
	cout << "We see the values are very close, but different." << endl;
	cout << "The absolute difference in the mean is: " << abs(mean(X_n) - mean(unif)) << endl;
	cout << "The absolute difference in the standard deviation is: " << abs(stdDev(X_n) - stdDev(unif)) << endl;
	cout << "We conclude the LGM algorithm performs very well." << endl;
	cout << "=========================================================================" << endl << endl;
	
	cout << "============================== Question 2  ==============================" << endl;
	cout << "============================== Question 2a ==============================" << endl;
	vector<int> bernTrials = indepBern(X_n, n);
	cout << "Computed the independent Bernoulli distribution" << endl;
	cout << "=========================================================================" << endl << endl;

	cout << "============================== Question 2b ==============================" << endl;
	cout << "The empirical mean of independent bernoulli: " << mean(bernTrials) << endl;
	cout << "The empirical standard deviation of independent bernoulli: " << stdDev(bernTrials) << endl;
	writeToTXT(bernTrials, "q2.txt");
	cout << "Successfully wrote to text q2 for R plotting" << endl;
	cout << "=========================================================================" << endl << endl;

	cout << "============================== Question 3  ==============================" << endl;
	cout << "============================== Question 3a ==============================" << endl;
	int bN = 44; // The binomial n
	double p = 0.64; // The binomial p
	int size = 1000;
	vector<double> binomial(size);
	binomial = Generator.rbin(size, bN, p);
	cout << "Generated the 1,000 random numbers of Binomial Distribution." << endl;
	cout << "Successfully wrote to text q3b for R plotting" << endl;
	cout << "=========================================================================" << endl << endl;
	cout << "============================== Question 3b ==============================" << endl;
	writeToTXT(binomial, "q3b.txt");
	int probabilityBin = probBin(binomial);
	cout << "Empirical estimates of P(X >= 40): " << probabilityBin / double(size) << endl;
	cout << "=========================================================================" << endl << endl;

	cout << "============================== Question 4  ==============================" << endl;
	cout << "============================== Question 4a ==============================" << endl;
	double lambda = 1.5;
	vector<double> exponential(n);
	exponential = Generator.rexp(n,lambda);
	cout << "Generated the 10,000 random numbers of the Exponential Distribution." << endl;
	writeToTXT(exponential, "q4c.txt");
	cout << "Successfully wrote to text q4c for R plotting" << endl;
	cout << "=========================================================================" << endl << endl;
	cout << "============================== Question 4b ==============================" << endl;
	// Counts the amount that is about 1,4
	int count1 = 0, count4 = 0;
	for (int i = 0; i < exponential.size(); i++) {
		if (exponential[i] >= 1) {
			count1 += 1;
		}
		if (exponential[i] >= 4) {
			count4 += 1;
		}
	}
	cout << "P(X >= 1): " << count1 / double(n) << endl;
	cout << "P(X >= 4): " << count4 / double(n) << endl;
	cout << "=========================================================================" << endl << endl;
	cout << "============================== Question 4c ==============================" << endl;
	cout << "The empirical mean of exponential: " << mean(exponential) << endl;
	cout << "The empirical standard deviation of exponential: " << stdDev(exponential) << endl;
	cout << "=========================================================================" << endl << endl;

	
	cout << "============================== Question 5  ==============================" << endl;
	cout << "============================== Question 5a ==============================" << endl;
	vector<double> boxMullerNorm(n);
	auto startTime1 = high_resolution_clock::now();
	boxMullerNorm = Generator.boxMuller(n);
	auto stopTime1 = high_resolution_clock::now();
	auto timeElapse1 = duration_cast<milliseconds> (stopTime1 - startTime1);
	cout << "The time (milliseconds) it took to run the Box Muller N(0,1): " << timeElapse1.count() << endl;
	cout << "=========================================================================" << endl << endl;

	cout << "============================== Question 5b ==============================" << endl;
	vector<double> polarMarsagliaNorm(n);
	auto startTime2 = high_resolution_clock::now();
	polarMarsagliaNorm = Generator.polarMarsaglia(n);
	auto stopTime2 = high_resolution_clock::now();
	auto timeElapse2 = duration_cast<milliseconds> (stopTime2 - startTime2);
	cout << "The time (milliseconds) it took to run the Polar Marsaglia N(0,1): " << timeElapse2.count() << endl;
	cout << "=========================================================================" << endl << endl;

	cout << "============================== Question 5c ==============================" << endl;
	n = 1000000;
	cout << "Using this many Normally distributed random number: " << n / 2.0 << endl;
	startTime1 = high_resolution_clock::now();
	boxMullerNorm = Generator.boxMuller(n);
	stopTime1 = high_resolution_clock::now();
	timeElapse1 = duration_cast<milliseconds> (stopTime1 - startTime1);
	cout << "The time (milliseconds) it took to run the Box Muller N(0,1): " << timeElapse1.count() << endl;
	
	startTime2 = high_resolution_clock::now();
	polarMarsagliaNorm = Generator.polarMarsaglia(n);
	stopTime2 = high_resolution_clock::now();
	timeElapse2 = duration_cast<milliseconds> (stopTime2 - startTime2);
	cout << "The time (milliseconds) it took to run the Polar Marsaglia N(0,1): " << timeElapse2.count() << endl;
	cout << "We see that the difference between the two methods is a few milliseconds in our initial run time" << endl;
	cout << "to be not siginificantly different. It is within computationally varied speed that we cannot compare." << endl;
	cout << "However, when the sequences get larger, the difference is speed will differ greatly." << endl;
	cout << "Polar Method is shown to be quicker method" << endl;
	cout << "=========================================================================" << endl << endl;
	return 0;
}

void printVector(vector<double> X) {
	for (int i = 0; i < X.size(); i++) {
		cout << X[i];
	}
	return;
}

template <class T>
double mean(vector<T> X) {
	double sum = 0;
	for (int i = 0; i < X.size(); i++) {
		sum += X[i];
	}
	sum /= X.size();
	return sum;
}

template <class T>
double stdDev(vector<T> X) {
	double mu = mean(X); // Mean of the sequence, mu
	double sum = 0;
	for (int i = 0; i < X.size(); i++) {
		sum += pow(X[i] - mu,2);
	}
	sum /= X.size();
	return sqrt(sum);
}

template <class T>
void writeToTXT(vector<T> arr, string filename) {
	ofstream file;
	file.open(filename);
	for (int i = 0; i < arr.size(); i++) {
		file << arr[i] << endl;
	}
	file.close();
}

vector<int> indepBern(vector<double> X, int n) {
	vector<int> y(n); // A vector to store the transformations
	for (int i = 0; i < X.size(); i++) {
		if(.3 >  X[i] && X[i] >= 0) {
			y[i] = -1;
		} 
		else if (.65  > X[i] && X[i] >= .3) {
			y[i] = 0;
		}
		else if (.85 > X[i] && X[i] >= .65) {
			y[i] = 1;
		}
		else if (1 >= X[i] && X[i] >=.85) {
			y[i] = 2;
		}
	}
	return y;
}


int probBin(vector<double> X) {
	int counter = 0;
	for (int i = 0; i < X.size(); i++) {
		if (X[i] >= 40) {
			counter++;
		}
	}
	return counter;
}