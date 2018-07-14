#include "NormalDistribution.h"

using namespace std;

static boost::random::mt19937 GLOBALGEN(49643);


NormalDistribution::NormalDistribution(double dMu, double dSigma )
	:m_dMu( 0.0 ), m_dSigma( 0.0 ), m_dLogNormFactor( 0.0 )
{
	SetParameters(dMu, dSigma);
}


NormalDistribution::~NormalDistribution()
{
}

/*Required by abstract base class; vParams[0] = mean, vParams[ 1 ] = Sigma*/
void NormalDistribution::SetParameters(const std::vector< double > & vParams)
{
	SetParameters(vParams[0], vParams[1]);

}

/* Set the paramters for the normal distirbution */
void NormalDistribution::SetParameters(double dMu, double dSigma)
{
	m_dMu = dMu;
	m_dSigma = dSigma;

	m_dLogNormFactor = log(1.0 / (sqrt(8.0*atan(1.0))*m_dSigma));

	//Use Boos wherever possible
	m_rngNormal = boost::random::normal_distribution<>(m_dMu, m_dSigma);	//Random number generation
	m_disNorm = boost::math::normal_distribution< >(m_dMu, m_dSigma);       //PDF, CDF, ect

}

vector< double > NormalDistribution::GetParameters() const
{
	vector<double> vParams(2, 0);
	GetParameters(vParams[0], vParams[1]);
	return(vParams);
}


void NormalDistribution::GetParameters(double &dMu, double &dSigma) const
{
	dMu = m_dMu;
	dSigma = m_dSigma;
}


double  NormalDistribution::PDF(double x) const
{
	return( boost::math::pdf(m_disNorm, x) );
}

double NormalDistribution::CDF(double x) const
{
	return( boost::math::cdf(m_disNorm, x) );
}

double  NormalDistribution::LogPDF(double x) const
{

	double dTmp = (x - m_dMu) / m_dSigma;
	dTmp *= -0.5*dTmp;
	return dTmp + m_dLogNormFactor;

}

/* Sample a value from the distribution */
double NormalDistribution::GetValue()  
{
	return m_rngNormal(GLOBALGEN);
}


