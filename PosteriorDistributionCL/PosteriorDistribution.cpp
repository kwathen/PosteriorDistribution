#include "PosteriorDistribution.h"
#include "RandomWalkMetropHast.h"
//using namespace MDACC_Biostat;
using namespace std;

using namespace boost::math;

PosteriorDistribution::PosteriorDistribution()
{
}


PosteriorDistribution::~PosteriorDistribution()
{
}

double PosteriorDistribution::operator()(const std::vector<double> & vParams)
{
	double dLogPriorPlusLogLike = CalcualteLogPrior(vParams) + CalculateLogLikelihood(vParams);
	return(dLogPriorPlusLogLike);
}

double PosteriorDistribution::CalcualteLogPrior(const std::vector<double> & vParams)
{

	//double dLogPrior = m_PriorMu.LogPDF(vParams[0]); // log(m_PriorMu.PDF(vParams[0]));
	double dLogPrior = log(pdf(m_PriorMu, vParams[0]));
	return(dLogPrior);

}

double PosteriorDistribution::CalculateLogLikelihood(const std::vector<double> & vParams)
{
	double dLogLikelihood = 0.0;


	for (int i = 0; i < m_nQtyDataPoints; ++i)
	{	
		dLogLikelihood +=  -0.5 * (m_vData[i] - vParams[0])  *  (m_vData[i] - vParams[0]) / vParams[1];
	}
	dLogLikelihood -= m_nQtyDataPoints * log(vParams[1]);
	return(dLogLikelihood);
}


void PosteriorDistribution::GetSample(int nQtySample, int nBurnIn, std::vector<double> &vSamples)
{
	RandomWalkMetropHast sampler(*this);

	sampler.Sample(nQtySample, nBurnIn, vSamples);

}

void PosteriorDistribution::SetData(std::vector<double> & vData)
{
	m_vData = vData;

	m_nQtyDataPoints = static_cast<int>( m_vData.size() );

}
void PosteriorDistribution::SetPrior(double dPriorMean, double dPriorSD)
{
	m_dPriorMean = dPriorMean;
	m_dPriorSD  = dPriorSD;

	m_PriorMu = normal(m_dPriorMean, m_dPriorSD);

}
