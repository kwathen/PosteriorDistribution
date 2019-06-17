#include "RandomWalkMetropHast.h"


#include "boost/random/normal_distribution.hpp"
#include "boost/random/mersenne_twister.hpp"

static boost::random::mt19937 GLOBALGEN(48209521);



using namespace std;
using namespace boost::math;

RandomWalkMetropHast::RandomWalkMetropHast( PosteriorDistribution & post )
{
	m_ptrPost = &post;
}


RandomWalkMetropHast::~RandomWalkMetropHast()
{
}


void RandomWalkMetropHast::Sample(int nQtySample, int nBurnIn, boost::numeric::ublas::matrix<double> &mSamples)
{


	//For Proof of concept just setting this to 2, which should come from the posterior distribution 
	int nQtyParams = 2;

	cout.precision(3);
	mSamples.resize(nQtySample, nQtyParams, false);
	SetJumpDistributions();


	boost::random::uniform_01<> unif;
	//Set the intial value

	boost::numeric::ublas::vector<double> vParams(nQtyParams, 0);
	boost::numeric::ublas::vector<double> vParamsProp(nQtyParams, 0);
	int i = 0;

	for (i = 0; i < nQtyParams; ++i)
	{
		vParams[i] = m_vJumpDist[i].GetValue();
	}


	//Make sure the Var > 0
	//vParams[1] = max(vParams[1], -vParams[1]);

	double dLogPriorPlusLogLike = m_ptrPost->CalculateLogPriorPlusLogLikelihood(vParams);
	double dLogPriorPlusLogLikeProp = 0.0;
	double dRatio = 0.0;

	int iParam = 0;
	vector< double > vAccept;

	vParamsProp = vParams; // Just to make sure the parameters are set.

	vAccept.resize(nQtyParams, 0.0); // Reset to 0
	for (i = 0; i < nBurnIn; ++i)
	{
		//One parameter at a time updating
		for (iParam = 0; iParam < nQtyParams; ++iParam)
		{
			vParamsProp[iParam] = vParams[iParam] + m_vJumpDist[iParam].GetValue();

			//TODO: Make this more general for checking bounds
			//if (iParam == 1)  //on the Std Dev param need to make sure > 0
			//	vParamsProp[1] = max(vParamsProp[1], 0.0001);

			dLogPriorPlusLogLikeProp = m_ptrPost->CalculateLogPriorPlusLogLikelihood(vParamsProp);  //Calculate the log-like + log-prior
			dRatio = exp(std::min(0.0, dLogPriorPlusLogLikeProp - dLogPriorPlusLogLike));
			dRatio = exp( dLogPriorPlusLogLikeProp - dLogPriorPlusLogLike);

			if (unif(GLOBALGEN) <= dRatio) //#Accept the move
			{
				vParams = vParamsProp;
				dLogPriorPlusLogLike = dLogPriorPlusLogLikeProp;
			}

		}

	}
	vAccept.resize(nQtyParams, 0.0); // Reset to 0
	for (i = 0; i < nQtySample; ++i)
	{
		//One parameter at a time updating
		for (iParam = 0; iParam < nQtyParams; ++iParam)
		{
			vParamsProp[iParam] = vParams[iParam] + m_vJumpDist[iParam].GetValue( );

			//TODO: Make this more general for checking bounds
			//if (iParam == 1)  //on the Std Dev param need to make sure > 0
			//	vParamsProp[1] = max(vParamsProp[1], 0.0001);

			dLogPriorPlusLogLikeProp = m_ptrPost->CalculateLogPriorPlusLogLikelihood(vParamsProp);  //Calculate the log-like + log-prior
			dRatio = exp(std::min(0.0, dLogPriorPlusLogLikeProp - dLogPriorPlusLogLike));
			dRatio = exp( dLogPriorPlusLogLikeProp - dLogPriorPlusLogLike);

			if (unif(GLOBALGEN) <= dRatio) //#Accept the move
			{
				vAccept[iParam]++;
				vParams = vParamsProp;
				dLogPriorPlusLogLike = dLogPriorPlusLogLikeProp;
			}

		}
		if( m_bPrint )
			cout << i << "\t" << vParams[0] << "\t" << vParams[1] << endl;
		
		//Save the params to be returned
		for (iParam = 0; iParam < nQtyParams; ++iParam)
		{
			mSamples(i, iParam) = vParams[iParam];
		}
		//vSamples[i] = vParams[0];

	}
	cout << "The acceptance rate for beta0 was " << vAccept[ 0 ]/ nQtySample << " and for beta1 " << vAccept[1] / nQtySample << endl;


}


void RandomWalkMetropHast::SetJumpDistributions()
{
	int i = 0;
	//For Proof of concept just setting this to 2, which should come from the posterior distribution 
	int nQtyParams = 2;
	m_vJumpDist.resize(2);
	for (i = 0; i < nQtyParams; ++i)
	{
		m_vJumpDist[i] = NormalDistribution(0, 2.0); //Starting jumping Std Dev is 2.0	
	}

	boost::random::uniform_01<> unif;
	//Set the intial value
	boost::numeric::ublas::vector<double> vParams(2, 0);
	boost::numeric::ublas::vector<double> vParamsProp(2, 0);

	for (i = 0; i < nQtyParams; ++i)
	{
		vParams[i] = m_vJumpDist[i].GetValue();
		vParamsProp[i] = vParams[i];
	}

	//Make sure the Var > 0
	//vParams[1] = max(vParams[1], -vParams[1]);
	vParamsProp[1] = vParams[1];
	double dLogPriorPlusLogLike = m_ptrPost->CalculateLogPriorPlusLogLikelihood(vParams);
	double dLogPriorPlusLogLikeProp = 0.0;

	//We are trying to set the jump std dev before we start the burn in and sampling of the chain
	//So we will interate through adjusting the jump up to 10 times and each time we will give a sample of 250 to see what the accpatace rate is
	//if all paramters have an acceptance rate > 0.2 and <0.6 then we are done.
	//TODO: Check the cutoffs of 0.2 and 0.6

	int nMaxAdj = 30;
	int nAdjChainLen = 50000;
	int iParam = 0;
	vector< double > vAccept;
	double dRatio = 0.0;
	
	bool bUpdateJump = false;
	vAccept.resize(nQtyParams, 0.0); // Reset to 0
	for (int j = 0; j < nMaxAdj; ++j)
	{
		vAccept[0] = 0;
		vAccept[1] = 0;
		for (i = 0; i < nAdjChainLen; ++i)
		{
			//One parameter at a time updating
			for (iParam = 0; iParam < nQtyParams; ++iParam)
			{
				vParamsProp[iParam] = vParams[iParam] + m_vJumpDist[iParam].GetValue( );
				//TODO: Make this more general for checking bounds
				//if (iParam == 1)  //on the Std Dev param need to make sure > 0
				//	vParamsProp[1] = max(vParamsProp[1], 0.0001);

				dLogPriorPlusLogLikeProp = m_ptrPost->CalculateLogPriorPlusLogLikelihood(vParamsProp);  //Calculate the log-like + log-prior
				dRatio = exp(std::min(0.0, dLogPriorPlusLogLikeProp - dLogPriorPlusLogLike));

				if (unif(GLOBALGEN) < dRatio) //#Accept the move
				{
					vParams = vParamsProp;
					dLogPriorPlusLogLike = dLogPriorPlusLogLikeProp;
					vAccept[iParam]++;
				}

			}

		}

		
		 bUpdateJump = false;
		//Check the acceptance and update the Std Dev of the jump variance if needed
		//TODO: This is not the optimal way, I had a more efficient way something like a bisection method to get to the desired accepatece rate

		double dSigma, dMu;
		for (iParam = 0; iParam < nQtyParams; ++iParam)
		{
			cout << "Parameter " << iParam << " Acceptance rate " << vAccept[iParam] / (1.0*nAdjChainLen) << "  ";

			if (vAccept[iParam] / (1.0*nAdjChainLen) < .15) // Need a smaller jump variance 
			{
				cout << "Decreasing Var" << endl;
				m_vJumpDist[iParam].GetParameters( dMu, dSigma );
				bUpdateJump = true;
				m_vJumpDist[iParam] =  NormalDistribution(dMu, dSigma /1.25);  //Cut the StdDev in half
			}
			else if (vAccept[iParam] / (1.0*nAdjChainLen) > 0.2) // Need a larger jump variance 
			{
				cout << "Increasing Var" << endl;
				m_vJumpDist[iParam].GetParameters(dMu, dSigma);
				bUpdateJump = true;
				m_vJumpDist[iParam] = NormalDistribution(dMu, dSigma * 1.15);  //Cut the StdDev in half
			}
		}
		if (bUpdateJump == false)
			break;
	}




}





