#include "PosteriorDistribution.h"
#include "RandomWalkMetropHast.h"

using namespace std;


//double PosteriorDistribution::operator()(const std::vector<double> & vParams)
//{
//	double dLogPriorPlusLogLike = CalcualteLogPrior(vParams) + CalculateLogLikelihood(vParams);
//	return(dLogPriorPlusLogLike);
//}
double PosteriorDistribution::CalculateLogPriorPlusLogLikelihood(const  boost::numeric::ublas::vector<double> & vParams)
{
	double dLogPriorPlusLogLike = CalcualteLogPrior(vParams) + CalculateLogLikelihood(vParams);
	return(dLogPriorPlusLogLike);
}

void PosteriorDistribution::GetSample(int nQtySample, int nBurnIn, boost::numeric::ublas::matrix<double> &mSamples)
{
	RandomWalkMetropHast sampler(*this);

	sampler.Sample(nQtySample, nBurnIn, mSamples);

}
