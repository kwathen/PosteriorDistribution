#pragma once

#include <boost/numeric/ublas/matrix.hpp> 
#include <boost/numeric/ublas/vector.hpp> 
#include "NormalDistribution.h"

class PosteriorDistribution
{
public:
	PosteriorDistribution() {};
	~PosteriorDistribution() {};
	//virtual double operator()(const std::vector<double>  & vParams);
	virtual double CalculateLogPriorPlusLogLikelihood(const boost::numeric::ublas::vector<double>  & vParams);
	virtual double CalcualteLogPrior(const boost::numeric::ublas::vector<double> & vParams) = 0;
	virtual double CalculateLogLikelihood(const boost::numeric::ublas::vector<double> & vParams) =0;


	virtual void SetData(struct & sData ) {};
	virtual void SetPrior(double dPriorMean, double dPriorSD) {};

	virtual void GetSample(int nQtySample, int nBurnIn, boost::numeric::ublas::matrix<double> &mSamples);

	

};
