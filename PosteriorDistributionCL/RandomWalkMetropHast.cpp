#include "RandomWalkMetropHast.h"


#include "boost/random/normal_distribution.hpp"
#include "boost/random/mersenne_twister.hpp"

static boost::random::mt19937 GLOBALGEN(48209521);



using namespace std;
using namespace boost::math;

RandomWalkMetropHast::RandomWalkMetropHast( PosteriorDistribution & post )
{
	m_Post = post;
}


RandomWalkMetropHast::~RandomWalkMetropHast()
{
}
//
//void RandomWalkMetropHast::Sample(int nQtySample, int nBurnIn, vector<double> &vSamples)
//{
//	vSamples.resize(nQtySample, 0.0);
//	boost::random::normal_distribution<> normMuProp(0, 0.6);
//	boost::random::normal_distribution<> normSigmaProp(0, 0.6);
//
//	boost::random::uniform_01<> unif;
//		//Set the intial value
//	double dMu = normMuProp(GLOBALGEN);
//
//	double dMuProp		= dMu;
//	double dSigma = max(normSigmaProp( GLOBALGEN ),0.0001) ;  //Can't be < 0 
//	double dSimgaProp	= dSigma;
//	double dRatio		= 0;
//
//	//Need to keep track of the current values of the parameters an the proposed parameters. 
//	vector<double> vParams(2, 0);
//	vector<double> vParamsProp(2, 0);
//	vParams[0] = dMu;
//	vParams[1] = dSigma;
//
//	
//	double dLogPriorPlusLogLike     = m_Post(vParams);
//	double dLogPriorPlusLogLikeProp = 0.0;
//
//	//Burn in 
//	int i = 0;
//	for (i = 0; i < nBurnIn; ++i)
//	{
//
//		vParamsProp[0] = vParams[0] + normMuProp(GLOBALGEN);
//		vParamsProp[1] = max(vParams[1] + normSigmaProp(GLOBALGEN), 0.01);
//		dLogPriorPlusLogLikeProp = m_Post(vParamsProp);
//
//		dRatio = exp(std::min(0.0, dLogPriorPlusLogLikeProp - dLogPriorPlusLogLike));
//
//		if (unif(GLOBALGEN) < dRatio ) //#Accept the move
//		{
//			vParams = vParamsProp;
//			dLogPriorPlusLogLike = dLogPriorPlusLogLikeProp;
//
//		}
//	}
//
//
//	//Now keep just vParam[0] eg the mean 
//	double dAcceptance1 = 0;
//	double dAcceptance2 = 0;
//	cout.precision(3);
//	for (i = 0; i < nQtySample; ++i)
//	{
//		//Propose both parameters
//		vParamsProp[0] = vParams[0] + normMuProp(GLOBALGEN);
//
//		dLogPriorPlusLogLikeProp = m_Post(vParamsProp);
//
//		dRatio = exp(std::min(0.0, dLogPriorPlusLogLikeProp - dLogPriorPlusLogLike));
//
//		cout << "mu " << vParams[0] << " " << vParamsProp[0] << " ratio of move " << dRatio << " ";
//		if (unif(GLOBALGEN) < dRatio ) //#Accept the move
//		{
//			vParams[0] = vParamsProp[0];
//			dLogPriorPlusLogLike = dLogPriorPlusLogLikeProp;
//			dAcceptance1++;
//
//		}
//		vSamples[i] = vParams[0];
//
//		vParamsProp[1] = max(vParams[1] + normSigmaProp(GLOBALGEN), 0.01);
//
//		dLogPriorPlusLogLikeProp = m_Post(vParamsProp);
//
//		dRatio = exp(std::min(0.0, dLogPriorPlusLogLikeProp - dLogPriorPlusLogLike));
//
//		cout << "\tsigma " << vParams[1] << " " << vParamsProp[1] << " ratio of move " << dRatio << " ";
//		if (unif(GLOBALGEN) < dRatio) //#Accept the move
//		{
//			vParams[1] = vParamsProp[1];
//			dLogPriorPlusLogLike = dLogPriorPlusLogLikeProp;
//			dAcceptance2++;
//
//		}
//		cout << i << " " << vParams[0] << " " << vParams[1] << endl;
//
//	}
//	cout << "The acceptance rate was " << dAcceptance1 / nQtySample << " " << dAcceptance2 / nQtySample << endl;
//	
//
//
//
//}


void RandomWalkMetropHast::Sample(int nQtySample, int nBurnIn, vector<double> &vSamples)
{


	cout.precision(3);

	vSamples.resize(nQtySample, 0.0);
	SetJumpDistributions();

	

	//For Proof of concept just setting this to 2, which should come from the posterior distribution 
	int nQtyParams = 2;

	boost::random::uniform_01<> unif;
	//Set the intial value

	vector<double> vParams(2, 0);
	vector<double> vParamsProp(2, 0);
	int i = 0;

	for (i = 0; i < nQtyParams; ++i)
	{
		vParams[i] = m_vJumpDist[i](GLOBALGEN);
	}


	//Make sure the Var > 0
	vParams[1] = max(vParams[1], -vParams[1]);

	double dLogPriorPlusLogLike = m_Post(vParams);
	double dLogPriorPlusLogLikeProp = 0.0;
	double dRatio = 0.0;

	int iParam = 0;
	vector< double > vAccept;

	vAccept.resize(nQtyParams, 0.0); // Reset to 0
	for (i = 0; i < nBurnIn; ++i)
	{
		//One parameter at a time updating
		for (iParam = 0; iParam < nQtyParams; ++iParam)
		{
			vParamsProp[iParam] = vParams[iParam] + m_vJumpDist[iParam](GLOBALGEN);

			//TODO: Make this more general for checking bounds
			if (iParam == 1)  //on the Std Dev param need to make sure > 0
				vParamsProp[1] = max(vParamsProp[1], 0.0001);

			dLogPriorPlusLogLikeProp = m_Post(vParamsProp);  //Calculate the log-like + log-prior
			dRatio = exp(std::min(0.0, dLogPriorPlusLogLikeProp - dLogPriorPlusLogLike));

			if (unif(GLOBALGEN) < dRatio) //#Accept the move
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
			vParamsProp[iParam] = vParams[iParam] + m_vJumpDist[iParam](GLOBALGEN);

			//TODO: Make this more general for checking bounds
			if (iParam == 1)  //on the Std Dev param need to make sure > 0
				vParamsProp[1] = max(vParamsProp[1], 0.0001);

			dLogPriorPlusLogLikeProp = m_Post(vParamsProp);  //Calculate the log-like + log-prior
			dRatio = exp(std::min(0.0, dLogPriorPlusLogLikeProp - dLogPriorPlusLogLike));

			if (unif(GLOBALGEN) < dRatio) //#Accept the move
			{
				vAccept[iParam]++;
				vParams = vParamsProp;
				dLogPriorPlusLogLike = dLogPriorPlusLogLikeProp;
			}

		}
		cout << i << "\t" << vParams[0] << "\t" << vParams[1] << endl;
		vSamples[i] = vParams[0];

	}
	cout << "The acceptance rate for mu was " << vAccept[ 0 ]/ nQtySample << " and for sigma " << vAccept[1] / nQtySample << endl;


}


void RandomWalkMetropHast::SetJumpDistributions()
{
	int i = 0;
	//For Proof of concept just setting this to 2, which should come from the posterior distribution 
	int nQtyParams = 2;
	m_vJumpDist.resize(2);
	for (i = 0; i < nQtyParams; ++i)
	{
		m_vJumpDist[i] = boost::random::normal_distribution<>(0, 2.0);  //Starting jumping Std Dev is 2.0
	}

	boost::random::uniform_01<> unif;
	//Set the intial value

	vector<double> vParams(2, 0);
	vector<double> vParamsProp(2, 0);

	for (i = 0; i < nQtyParams; ++i)
	{
			vParams[i] = m_vJumpDist[i](GLOBALGEN);
			vParamsProp[i] = vParams[i];
	}

	//Make sure the Var > 0
	vParams[1] = max(vParams[1], -vParams[1]);
	vParamsProp[1] = vParams[1];
	double dLogPriorPlusLogLike = m_Post(vParams);
	double dLogPriorPlusLogLikeProp = 0.0;

	//We are trying to set the jump std dev before we start the burn in and sampling of the chain
	//So we will interate through adjusting the jump up to 10 times and each time we will give a sample of 250 to see what the accpatace rate is
	//if all paramters have an acceptance rate > 0.2 and <0.6 then we are done.
	//TODO: Check the cutoffs of 0.2 and 0.6

	int nMaxAdj = 10;
	int nAdjChainLen = 250;
	int iParam = 0;
	vector< double > vAccept;
	double dRatio = 0.0;
	
	bool bUpdateJump = false;
	for (int j = 0; j < nMaxAdj; ++j)
	{
		vAccept.resize(nQtyParams, 0.0); // Reset to 0
		for (i = 0; i < nAdjChainLen; ++i)
		{
			//One parameter at a time updating
			for (iParam = 0; iParam < nQtyParams; ++iParam)
			{
				vParamsProp[iParam] = vParams[iParam] + m_vJumpDist[iParam](GLOBALGEN);

				//TODO: Make this more general for checking bounds
				if (iParam == 1)  //on the Std Dev param need to make sure > 0
					vParamsProp[1] = max(vParamsProp[1], 0.0001);

				dLogPriorPlusLogLikeProp = m_Post(vParamsProp);  //Calculate the log-like + log-prior
				dRatio = exp(std::min(0.0, dLogPriorPlusLogLikeProp - dLogPriorPlusLogLike));

				if (unif(GLOBALGEN) < dRatio) //#Accept the move
				{
					vParams = vParamsProp;
					dLogPriorPlusLogLike = dLogPriorPlusLogLikeProp;
					vAccept[iParam]++;
				}

			}

		}


		bool bUpdateJump = false;
		//Check the acceptance and update the Std Dev of the jump variance if needed
		//TODO: This is not the optimal way, I had a more efficient way something like a bisection method to get to the desired accepatece rate

		for (iParam = 0; iParam < nQtyParams; ++iParam)
		{
			
			if (vAccept[iParam] / nAdjChainLen < 0.2) // Need a smaller jump variance 
			{
				bUpdateJump = true;
				m_vJumpDist[iParam] =  boost::random::normal_distribution<>(0, m_vJumpDist[iParam].sigma() /2.0);  //Cut the StdDev in half
			}
			else if (vAccept[iParam] / nAdjChainLen > 0.6) // Need a larger jump variance 
			{
				bUpdateJump = true;
				m_vJumpDist[iParam] = boost::random::normal_distribution<>(0, m_vJumpDist[iParam].sigma() * 1.75);  //Cut the StdDev in half
			}
		}
		if (bUpdateJump == false)
			break;
	}




}





