#pragma once

#include "PosteriorDistribution.h"

#include <boost/numeric/ublas/matrix.hpp> 
#include <boost/numeric/ublas/vector.hpp> 

class LogisticRegressionPosterior :public PosteriorDistribution
{
public:
	LogisticRegressionPosterior();
	~LogisticRegressionPosterior();

	double CalcualteLogPrior(const boost::numeric::ublas::vector<double> & vParams);
	double CalculateLogLikelihood(const boost::numeric::ublas::vector<double> & vParams);

	void SetData( const std::vector<double> & vDose, const std::vector<int> vEvent );
	void SetPrior(const boost::numeric::ublas::matrix<double>& mPriorMean, 
		          const boost::numeric::ublas::matrix<double> & mVarCov );

private:
	std::vector<double>			m_vDose;
	std::vector< int >			m_vEvent;
	int                         m_nQtyDataPoints;
	boost::numeric::ublas::matrix<double>  m_mPriorMean;
	boost::numeric::ublas::matrix<double>  m_mPriorVarCov;
	boost::numeric::ublas::matrix<double>  m_mPriorVarCovInv;

};

