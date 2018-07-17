#include "Functions.h"

bool InvertMatrix(const boost::numeric::ublas::matrix<double>& input, boost::numeric::ublas::matrix<double>& inverse)
{
	using namespace boost::numeric::ublas;
	typedef permutation_matrix<std::size_t> pmatrix;

	matrix<double> A(input);
	pmatrix pm(A.size1());
	int res = lu_factorize(A, pm);

	if (res != 0) 
		return false;
	inverse.resize(input.size1(), input.size2());
	inverse.assign(identity_matrix<double>(A.size1()));

	lu_substitute(A, pm, inverse);
	return true;
}

//value, mean vector and the var-cov matrix
double MultivariateNormalLogPDF(const boost::numeric::ublas::vector<double>& vX,
	                            const boost::numeric::ublas::matrix<double>& mMean,
	                            const boost::numeric::ublas::matrix<double> & mVarCov,
								bool bIsInverse )
{
	namespace ublas = boost::numeric::ublas;
	using namespace ublas;
	boost::numeric::ublas::matrix<double> mX(vX.size(), 1);
	column(mX, 0) = vX;
	

	boost::numeric::ublas::matrix<double> mXMinusMean;
	mXMinusMean = mX - mMean;

	//ublas::matrix<double> vec( dim, 1);
	//ublas::matrix<double> mat(dim, dim);
	ublas::matrix<double> mVarCovInv( mVarCov.size1(), mVarCov.size2() );
	if (bIsInverse)
	{
		mVarCovInv = mVarCov;  //The matrix that was supplied was the inverse, no need to compute.  
	}
	else
	{
		InvertMatrix(mVarCov, mVarCovInv);
	}
	
	boost::numeric::ublas::matrix<double>  vec1 = prod( trans(mXMinusMean), mVarCovInv);
	double dLogPDF = -0.5 * prod(vec1, mXMinusMean)(0,0);
	return(dLogPDF);


}

