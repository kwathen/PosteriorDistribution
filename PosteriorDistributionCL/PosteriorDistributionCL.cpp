// PosteriorDistributionCL.cpp : Defines the entry point for the console application.
//
#include <boost/lambda/lambda.hpp>
#include "boost/math/distributions/normal.hpp"
#include <iostream>
#include <iterator>
#include <numeric>
#include <algorithm>
#include <cmath>
#include<fstream>
#include <iostream>


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

#include "SimpleNormalPosterior.h"
#include "LogisticRegressionPosterior.h"
#include "NormalDistribution.h"

#include "Functions.h"
using namespace std;


using namespace boost::math;


//static boost::random::mt19937 GLOBALGEN(48209521);
static boost::random::mt19937 GLOBALGEN;// (443);



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

int main()
{
	//Simulate the example data - using the NormalDistribtuon class
	//int nSampSize = 200;
	//NormalDistribution normal(10, 3);
	//std::vector<double> vData = normal.GetValues(nSampSize);
	////boost::numeric::ublas::vector<double> vData


	////For a Dll the following code probably is what we would want to wrap and make callable from R

	//boost::numeric::ublas::matrix<double> mSamples;
	//SimpleNormalPosterior post;
	//post.SetData(vData);
	//post.SetPrior(0, 100);

	////When sampling This will now print the sampled values as  mu, sigma 
	//post.GetSample(1000, 1000, mSamples);


	//cout << "Sample mean of the observed data = " << SampleMean(vData) << endl;
	//cout << "Posterior mean = " << SampleMean( column( mSamples, 0 ) ) << endl;


	//// An example of things can can be done with the NormalDistribution
	//cout << endl << " Example using the NormalDistribution Class ..." << endl;
	//cout << "PDF(10 ) = " << normal.PDF(10) << " CDF(10) = " << normal.CDF(10) << " LogPFD( 10 ) = " << normal.LogPDF(10);
	//cout << " Example Value " << normal.GetValue() << endl;

	//Setup for the BLRM
	boost::numeric::ublas::matrix<double> mSamples;
	boost::numeric::ublas::matrix<double> mVarCov(2, 2);
	boost::numeric::ublas::matrix<double> mMean(2, 1);
	boost::numeric::ublas::vector<double> vBeta(2);

	//Set the prior
	mVarCov(0, 0) = 4.0;
	mVarCov(1, 1) = 1.0;
	mVarCov(0, 1) = mVarCov(1, 0) = 0.0;

	mMean(0, 0) = -2.944439;
	mMean(1, 0) = 0;

	//Example betas to calcuate the log prior and log likelihood
	vBeta(0) = -1.2;
	vBeta(1) = -0.5;

	//Example data - same as in R 
	std::vector<double> vDose(3, 0);
	std::vector< int > vEvent(3, 0);

	vDose[0] = vDose[1] = vDose[2] = log( 0.3/2.4);
	//vDose[3] = vDose[4] = vDose[5] = log(0.6 / 2.4);
	/*
	vDose[3] = vDose[4] = vDose[5] = 1.0;
	vDose[6] = vDose[7] = vDose[8] = 2.0;

	vDose[0] = vDose[1] = vDose[2] = 0.0;
	vDose[3] = vDose[4] = vDose[5] = 1.0;
	vDose[6] = vDose[7] = vDose[8] = 2.0;

	vEvent[7] = vEvent[8] = 1;
	vEvent[3] = 1;
	
	*///Sample the posterior
	vEvent[2] = 1;
	//vEvent[4] = 1;
 	LogisticRegressionPosterior blrmPost;
	
	blrmPost.SetData(vDose, vEvent);
	blrmPost.SetPrior(mMean, mVarCov);

	cout << "Compare to the computed values in R" << endl;
	cout << "Log Prior = " << blrmPost.CalcualteLogPrior(vBeta) << endl;
	cout << "Log Like = " << blrmPost.CalculateLogLikelihood(vBeta) << endl;

	double dPrGr0 = 0.0;
	double dMean0 = 0.0;
	double dMean1 = 0.0;

	cout << "Chain 1 ***" << endl;
	int nQtySampe = 2000000;
//	blrmPost.GetSample(nQtySampe, 1000000, mSamples);
	
	std::ofstream outfile;
	outfile.open("C:\\Kyle\\OneDrive\\Contract Work\\AlexD\\PosteriorDistribution\\Ch2Pack\\sample.txt");
	outfile.precision(4);
	outfile << "msamps1 <- matrix(c(";

	dPrGr0 = 0;
/*
	for (int i = 0; i < nQtySampe; i+=20)
	{
		outfile << mSamples(i, 0) << ", " << mSamples(i, 1) << ", ";
		if (mSamples(i, 1) > 0.0)
			dPrGr0++;
		dMean0 += mSamples(i, 0);
		dMean1 += mSamples(i, 1);
	}
	cout  << "Pr( beta1 > 0 |data ) = " << dPrGr0 / nQtySampe << " mean0 " << dMean0/ nQtySampe << " mean1 " << dMean1/ nQtySampe << endl << endl;


	cout << "Chain 2 ***" << endl;
	blrmPost.GetSample(nQtySampe, 100000, mSamples);

	dPrGr0 =  dMean0 = dMean1 = 0.0;
	for (int i = 0; i < nQtySampe; i += 20)
	{
		outfile << mSamples(i, 0) << ", " << mSamples(i, 1) << ", ";
		if (mSamples(i, 1) > 0.0)
			dPrGr0++;
		dMean0 += mSamples(i, 0);
		dMean1 += mSamples(i, 1);

	}
	cout << "Pr( beta1 > 0 |data ) = " << dPrGr0 / nQtySampe << " mean0 " << dMean0 / nQtySampe << " mean1 " << dMean1 / nQtySampe << endl << endl;
*/


	cout << "Chain 3 ***" << endl;
	blrmPost.GetSample(10*nQtySampe, 300000, mSamples);
	

	dPrGr0 = dMean0 = dMean1 = 0.0;
	for (int i = 0; i < nQtySampe; i += 10)
	{
		if( i + 10 < nQtySampe)
			outfile << mSamples(i, 0) << ", " << mSamples(i, 1) << ", ";
		else
			outfile << mSamples(i, 0) << ", " << mSamples(i, 1) << "), byrow=TRUE, ncol=2) ";
		if (mSamples(i, 1) > 0.0)
			dPrGr0++;
		dMean0 += mSamples(i, 0);
		dMean1 += mSamples(i, 1);

	}
	outfile.close();
	cout << "Pr( beta1 > 0 |data ) = " << dPrGr0 / nQtySampe << " mean0 " << dMean0 / nQtySampe << " mean1 " << dMean1 / nQtySampe << endl << endl;

	return(0);

}
