#pragma once

//#include <BiostatGeneral.h>
#include <vector>

#include "boost/math/distributions/normal.hpp"


class PosteriorDistribution
{
public:
	PosteriorDistribution();
	~PosteriorDistribution();
	double operator()(const std::vector<double>  & vParams);
	double CalcualteLogPrior(const std::vector<double> & vParams);
	double CalculateLogLikelihood(const std::vector<double> & vParams);

	void SetData(std::vector<double> & vData);
	void SetPrior(double dPriorMean, double dPriorSD);
	void GetSample(int nQtySample, int nBurnIn, std::vector<double> &vSamples);

private:
	std::vector<double>			m_vData;
	double						m_dPriorMean;
	double						m_dPriorSD;
	int                         m_nQtyDataPoints;
	boost::math::normal_distribution<> m_PriorMu;

};
