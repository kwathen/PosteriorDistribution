#include "LogisticRegressionPosterior.h"
#include "Functions.h"

#include <iostream>
#include <cmath>

using namespace std;
LogisticRegressionPosterior::LogisticRegressionPosterior()
{
}


LogisticRegressionPosterior::~LogisticRegressionPosterior()
{
}

double LogisticRegressionPosterior::CalcualteLogPrior(const boost::numeric::ublas::vector<double> & vParams)
{
	return(MultivariateNormalLogPDF(vParams, m_mPriorMean, m_mPriorVarCovInv, true));
}

double LogisticRegressionPosterior::CalculateLogLikelihood(const boost::numeric::ublas::vector<double> & vParams)
{
	double dLogLike = 0.0;
	double dEta = 0.0;
	for (int i = 0; i < m_nQtyDataPoints; ++i)
	{
		dEta = vParams[0] + m_vDose[i] * exp(vParams[1]);

		if (m_vEvent[i] == 1)  //Outcome is a 1, add the log( p )
		{
			if (dEta < 0)
			{
				dLogLike += (dEta - log1p(exp(dEta)));
			}
			else
			{
				dLogLike -= log1p(exp(-dEta));
			}

		}
		else if(m_vEvent[i] == 0 )
		{
			if (dEta < 0)
			{
				dLogLike -=  log1p(exp(dEta));
			}
			else
			{
				dLogLike -= ( dEta + log1p( exp( -dEta ) ) );
			}

		}
		else
		{
			//This is an error -this should stop execution
			cout << "ERROR IN THE DATA********************" << endl;
		}

	}
	return(dLogLike);
}



void LogisticRegressionPosterior::SetData(const std::vector<double> & vDose, const std::vector<int> vEvent)
{
	m_vDose = vDose;
	m_vEvent = vEvent;
	m_nQtyDataPoints = static_cast<int>( vDose.size() );
	

}
void LogisticRegressionPosterior::SetPrior(	const boost::numeric::ublas::matrix<double>& mPriorMean,
											const boost::numeric::ublas::matrix<double> & mPriorVarCov)
{
	
	m_mPriorMean = mPriorMean;
	m_mPriorVarCov = mPriorVarCov;
	InvertMatrix(mPriorVarCov, m_mPriorVarCovInv );
}