#pragma once

#include <vector>
#include "NormalDistribution.h"
#include "PosteriorDistribution.h"
//#include "boost/random/normal_distribution.hpp"

class RandomWalkMetropHast
{
public:
	RandomWalkMetropHast():m_bPrint(false) {}
	RandomWalkMetropHast( PosteriorDistribution & post );
	~RandomWalkMetropHast();
	void Sample(int nQtySample, int nBurnIn, boost::numeric::ublas::matrix<double> &mSamples);

	bool m_bPrint;
private:

	void SetJumpDistributions( );
	std::vector< NormalDistribution > m_vJumpDist;
	//std::vector< boost::random::normal_distribution<> > m_vJumpDist;

	PosteriorDistribution * m_ptrPost;
	//boost::math::normal m_Proposal;
	NormalDistribution m_Proposal;

	

};



