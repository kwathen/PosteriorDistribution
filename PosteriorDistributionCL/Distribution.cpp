#include "Distribution.h"


using namespace std;
Distribution::Distribution()
{
}


Distribution::~Distribution()
{
}

double Distribution::LogPDF(double x) const 
{
	return(log(PDF(x)));
}


std::vector< double >  Distribution::GetValues(int nQtySample) 
{
	vector<double> vValues(nQtySample, 0.0);
	for (int i = 0; i < nQtySample; ++i)
	{
		vValues[i] = GetValue( );
	}
	return( vValues );
}