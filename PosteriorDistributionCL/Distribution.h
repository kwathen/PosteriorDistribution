#pragma once

#include <vector>


class Distribution
{
public:
	Distribution();
	~Distribution();

	virtual void SetParameters(const std::vector< double > & vParams) = 0;  //Pure virtual so you must define for any derived class
	
	virtual std::vector< double > GetParameters() const = 0;

	virtual double PDF(double x) const = 0;                                 // const = 0 --> Pure virtual so you must define, const --> won't change distribution
	virtual double CDF(double x) const = 0;                                 // const = 0 --> Pure virtual so you must define, const --> won't change distribution

	/*! 
	You could use LogPDF(x) = log(PDF(x)) but there since many distributions have exp(f(x)) 
	it is better to just compute f(x) rather than log( exp( f(x) )) which is wasteful, and may fail due to underflow.
	*/
	virtual double LogPDF(double x) const ;


	/* Sample a value from the distribution */
	virtual double GetValue() { return 0;}


	/* Sample a vector of values from the distribution */
	virtual std::vector< double >  GetValues(int nQtySample);


};

