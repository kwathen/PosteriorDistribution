#pragma once
#include <boost/numeric/ublas/matrix.hpp> 
#include <boost/numeric/ublas/vector.hpp> 
#include <boost/numeric/ublas/io.hpp> 
#include <boost/numeric/ublas/vector_proxy.hpp> 
#include <boost/numeric/ublas/triangular.hpp> 
#include <boost/numeric/ublas/lu.hpp> 

//value, mean vector and the var-cov matrix

bool InvertMatrix(const boost::numeric::ublas::matrix<double>& input, boost::numeric::ublas::matrix<double>& inverse);

double MultivariateNormalLogPDF(const boost::numeric::ublas::vector<double>& vX,
	const boost::numeric::ublas::matrix<double>& mMean,
	const boost::numeric::ublas::matrix<double> & mVarCov,
	bool bIsInverse = false);