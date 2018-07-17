#include "SimpleNormalPosterior.h"

#include "RandomWalkMetropHast.h"

using namespace std;

SimpleNormalPosterior::SimpleNormalPosterior()
{
}


SimpleNormalPosterior::~SimpleNormalPosterior()
{
}



double SimpleNormalPosterior::CalcualteLogPrior(const boost::numeric::ublas::vector<double> & vParams)
{

	double dLogPrior = m_PriorMu.LogPDF(vParams[0]);
	return(dLogPrior);

}

//This function is to return the log( likelihood )
double SimpleNormalPosterior::CalculateLogLikelihood(const boost::numeric::ublas::vector<double> & vParams)
{
	double dLogLikelihood = 0.0;

	//There are two options for doing this -- which is easier/clear/less likely to have bugs?

	//Option 1 - This approach uses the NormalDistribution and is generic to other models,
	//           however, this approach is less efficient.  
	NormalDistribution likelihood(vParams[0], vParams[1]);
	for (int i = 0; i < m_nQtyDataPoints; ++i)
	{
		dLogLikelihood += likelihood.LogPDF(m_vData[i]);
	}


	//Option 2 - This approach uses is more efficient, but less clear.  I believe doing this is more
	//           likeley to have bugs as it did when I wrote it the first time and forgot to square 
	//for (int i = 0; i < m_nQtyDataPoints; ++i)
	//{	
	//	dLogLikelihood +=  -0.5 * (m_vData[i] - vParams[0])  *  (m_vData[i] - vParams[0]) / vParams[1];
	//}
	//dLogLikelihood -= m_nQtyDataPoints * log(vParams[1]);

	return(dLogLikelihood);
}


void SimpleNormalPosterior::SetData(std::vector<double> & vData)
{
	m_vData = vData;
	m_nQtyDataPoints = static_cast<int>(m_vData.size());

}
void SimpleNormalPosterior::SetPrior(double dPriorMean, double dPriorSD)
{
	m_dPriorMean = dPriorMean;
	m_dPriorSD = dPriorSD;


	m_PriorMu.SetParameters(m_dPriorMean, m_dPriorSD);

}
