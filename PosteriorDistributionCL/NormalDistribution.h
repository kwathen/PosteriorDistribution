#pragma once
#include "Distribution.h"

#include "boost/math/distributions/normal.hpp"
#include "boost/random/normal_distribution.hpp"
#include "boost/random/mersenne_twister.hpp"


class NormalDistribution :public Distribution
{
	public:
		NormalDistribution( double dMu = 0.0, double dSigma = 1.0 );
		~NormalDistribution();


		/*Required by abstract base class; vParams[0] = mean, vParams[ 1 ] = Sigma*/
		virtual void SetParameters(const std::vector< double > & vParams) ; 

		/* Set the paramters for the normal distirbution */
		virtual void SetParameters(double dMu, double dSigma);

		virtual std::vector< double > GetParameters() const;
		virtual void GetParameters(double &dMu, double &dSigma) const;

		virtual double PDF(double x) const;                                 // const won't change distribution
		virtual double CDF(double x) const;                                 // const won't change distribution


		/*!
		You could use LogPDF(x) = log(PDF(x)) but there since many distributions have exp(f(x))
		it is better to just compute f(x) rather than log( exp( f(x) )) which is wasteful, and may fail due to underflow.
		*/
		virtual double LogPDF(double x) const;


		/* Sample a value from the distribution */
		virtual double GetValue();


	private:
		double		m_dMu;				// Mean
		double		m_dSigma;			// Stadard Deviation
		double      m_dLogNormFactor;	//Log normal factor used for computing the LogPDF

		//Using Boost for random number generation, PDF ect
		boost::random::normal_distribution<> m_rngNormal;
		boost::math::normal_distribution<> m_disNorm;
};

