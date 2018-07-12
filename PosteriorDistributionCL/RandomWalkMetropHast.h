#pragma once

#include <vector>

#include "PosteriorDistribution.h"
#include "boost/random/normal_distribution.hpp"

class RandomWalkMetropHast
{
public:
	RandomWalkMetropHast( PosteriorDistribution & post );
	~RandomWalkMetropHast();
	void Sample(int nQtySample, int nBurnIn, std::vector<double> &vSamples);


private:

	void SetJumpDistributions( );
	std::vector< boost::random::normal_distribution<> > m_vJumpDist;

	PosteriorDistribution m_Post;
	boost::math::normal m_Proposal;

	

};


