// PosteriorDistributionCL.cpp : Defines the entry point for the console application.
//
#include <boost/lambda/lambda.hpp>
#include "boost/math/distributions/normal.hpp"
#include <iostream>
#include <iterator>
#include <numeric>
#include <algorithm>


#include "boost/math/distributions/students_t.hpp"
#include "boost/math/distributions/normal.hpp"
#include "boost/math/distributions/chi_squared.hpp"
#include "boost/math/distributions/beta.hpp"
#include "boost/math/distributions/gamma.hpp"

#include "boost/random/mersenne_twister.hpp"
#include "boost/random/uniform_real_distribution.hpp"
#include "boost/random/exponential_distribution.hpp"
#include "boost/random/normal_distribution.hpp"
#include "boost/random/beta_distribution.hpp"
#include "boost/math/distributions/normal.hpp"

#include "PosteriorDistribution.h"
#include "NormalDistribution.h"

using namespace std;


//static boost::random::mt19937 GLOBALGEN(48209521);
static boost::random::mt19937 GLOBALGEN(4289643);


// Vector of normally distributed values 
vector<double> Normal(const int &n, const double &mean, const double &sd) {

	vector<double> result(n);
	boost::random::normal_distribution<> distribution(mean, sd);
	for (int i = 0; i < n; i++) {
		result[i] = distribution(GLOBALGEN);
	}

	return result;

}

// Vector of uniformly distributed values 
vector<double> Uniform(const int &n, const double &min, const double &max) {

	std::vector<double> result(n);
	boost::random::uniform_real_distribution<> distribution(min, max);
	for (int i = 0; i < n; i++) {
		result[i] = distribution(GLOBALGEN);
	}

	return result;

}

//Function to compute the mean of a vector
double SampleMean( const  std::vector<double> & vObs)
{
	return(std::accumulate(vObs.begin(), vObs.end(), 0) / ( 1.0* vObs.size() ));
}


int main()
{
	//Simulate the example data - using the NormalDistribtuon class
	int nSampSize = 200;
	NormalDistribution normal(10, 3);
	vector<double> vData = normal.GetValues(nSampSize);


	//For a Dll the following code probably is what we would want to wrap and make callable from R

	vector<double> vSamples;
	PosteriorDistribution post;
	post.SetData(vData);
	post.SetPrior(0, 100);

	//When sampling This will now print the sampled values as  mu, sigma 
	post.GetSample(1000, 1000, vSamples);


	cout << "Sample mean of the observed data = " << SampleMean(vData) << endl;
	cout << "Posterior mean = " << SampleMean( vSamples ) << endl;


	// An example of things can can be done with the NormalDistribution
	cout << endl << " Example using the NormalDistribution Class ..." << endl;
	cout << "PDF(10 ) = " << normal.PDF(10) << " CDF(10) = " << normal.CDF(10) << " LogPFD( 10 ) = " << normal.LogPDF(10);
	cout << " Example Value " << normal.GetValue() << endl;
    return 0;
}
