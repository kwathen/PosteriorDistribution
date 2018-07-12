// PosteriorDistributionCL.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <boost/lambda/lambda.hpp>
#include "boost/math/distributions/normal.hpp"
#include <iostream>
#include <iterator>
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

using namespace std;


//static boost::random::mt19937 GLOBALGEN(48209521);
static boost::random::mt19937 GLOBALGEN(49643);


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

#include "boost/math/distributions/normal.hpp"
#include "PosteriorDistribution.h"
int main()
{

	int nSampSize = 200;
	vector<double> vData = Normal(nSampSize, 2.0, 2.0);


	vector<double> vSamples;
	PosteriorDistribution post;
	post.SetData(vData);
	post.SetPrior(0, 100);

	//When sampling This will now print the sampled values as  mu, sigma 
	post.GetSample(1000, 1000, vSamples);


	double dSum = 0.0;
	for (int i = 0; i < nSampSize; ++i)
	{
		//cout << vData[i] << endl;
		dSum += vData[i];
	}
	cout << "Sample mean of the observed data = " << dSum / nSampSize << endl;

	dSum = 0.0;
	for (int i = 0; i < 1000; ++i)
	{
		dSum += vSamples[i];
	}

	cout << "Posterior mean = " << dSum / 1000.0 << endl;

    return 0;
}

/* 
Example code from boost
double dval = 2.4;
boost::math::normal_distribution<> distribution(2, 3);

double dPDF = boost::math::pdf(distribution, dval);




cout << "Unfiorms " << endl;
boost::random::uniform_01<> unif;
cout << unif(GLOBALGEN) << " " << unif(GLOBALGEN) << endl;


*/