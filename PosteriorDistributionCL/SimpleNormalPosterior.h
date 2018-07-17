#pragma once
#include "PosteriorDistribution.h"

#include <vector>
#include "NormalDistribution.h"

class SimpleNormalPosterior : public PosteriorDistribution
{
public:
	SimpleNormalPosterior();
	~SimpleNormalPosterior();

	double CalcualteLogPrior(const boost::numeric::ublas::vector<double> & vParams);
	double CalculateLogLikelihood(const boost::numeric::ublas::vector<double> & vParams);

	void SetData(std::vector<double> & vData);
	void SetPrior(double dPriorMean, double dPriorSD);

private:
	std::vector<double>			m_vData;
	double						m_dPriorMean;
	double						m_dPriorSD;
	int                         m_nQtyDataPoints;
	NormalDistribution          m_PriorMu;
};
