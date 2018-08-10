#include <cassert>
#include <limits.h>
#include <complex>
#include <unsupported/Eigen/MatrixFunctions>
#include "matrix.h"

using namespace linopt;

complex_type linopt::permanent(const matrix_type &M)
{
	assert(M.cols() == M.rows());
	if(M.cols() == 2)
		return M(0, 0)*M(1, 1) + M(0, 1)*M(1, 0);
	if(M.cols() == 3)
		return M(0, 0)*(M(1, 1)*M(2, 2) + M(1, 2)*M(2, 1)) +
			   M(0, 1)*(M(1, 0)*M(2, 2) + M(1, 2)*M(2, 0)) +
			   M(0, 2)*(M(1, 0)*M(2, 1) + M(1, 1)*M(2, 0));
	if(M.cols() == 4)
		return M(0,1)*M(1,3)*M(2,2)*M(3,0) + M(0,1)*M(1,2)*M(2,3)*M(3,0) +
			   M(0,0)*M(1,3)*M(2,2)*M(3,1) + M(0,0)*M(1,2)*M(2,3)*M(3,1) +
			   M(0,1)*M(1,3)*M(2,0)*M(3,2) + M(0,0)*M(1,3)*M(2,1)*M(3,2) +
			   M(0,1)*M(1,0)*M(2,3)*M(3,2) + M(0,0)*M(1,1)*M(2,3)*M(3,2) +
			   M(0,3)*(M(1,2)*(M(2,1)*M(3,0) + M(2,0)*M(3,1)) +
			   M(1,1)*(M(2,2)*M(3,0) + M(2,0)*M(3,2)) + M(1,0)*(M(2,2)*M(3,1) +
			   M(2,1)*M(3,2))) + M(0,1)*M(1,2)*M(2,0)*M(3,3) +
			   M(0,0)*M(1,2)*M(2,1)*M(3,3) + M(0,1)*M(1,0)*M(2,2)*M(3,3) +
			   M(0,0)*M(1,1)*M(2,2)*M(3,3) + M(0,2)*(M(1,3)*(M(2,1)*M(3,0) +
			   M(2,0)*M(3,1)) + M(1,1)*(M(2,3)*M(3,0) + M(2,0)*M(3,3)) +
			   M(1,0)*(M(2,3)*M(3,1) + M(2,1)*M(3,3)));
	assert(static_cast<unsigned long>(M.cols()) <= CHAR_BIT*sizeof(unsigned long));
	complex_type perm = 0.;
	vector_type sum = M.rowwise().sum();
	real_type mult;
	unsigned long n = 0ul;
	unsigned long nmax = 1ul << (M.cols() - 1ul);
	while(n < nmax)
	{
		complex_type sp = sum.prod();
		perm += (n & 1ul) ? -sp : sp;
		unsigned long i = n ^ (n + 1ul);
		i++;
		n++;
		mult = (n & i) ? 2 : -2;
		i = __builtin_ctz(i) - 1ul;
		sum += mult * M.col(i);
	}
	perm /= nmax;
	return perm;
}

void linopt::hurwitz_parametrization(matrix_type &M, const point &x)
{
	const int N = std::sqrt(x.size());
	assert(N*N == static_cast<int>(x.size()));
	if(M.cols() != N || M.rows() != N)
		M.resize(N, N);
	int i, j, k = 0;
	complex_type eii, eij;
	vector_type coli(N);
	real_type xi;
	eii = std::polar(1., 2.*M_PI*x[k++]);
	M.setZero();
	M.diagonal().fill(eii);
	for(j = 1; j <= N-1; j++)
	{
		for(i = j - 1; i >= 1; i--)
		{
			xi = pyramid(x[k++], 1.);
			xi = std::pow(xi, 1./(2.*(i + 1)));
			eii = std::polar(std::sqrt(1. - xi*xi), 2.*M_PI*x[k++]);
			eij = xi;
			coli = M.col(i);
			M.col(i) = M.col(i)*eii - M.col(j)*conj(eij);
			M.col(j) =     coli*eij + M.col(j)*conj(eii);
		}
		xi = pyramid(x[k++], 1.);
		xi = std::sqrt(xi);
		eii = std::polar(std::sqrt(1. - xi*xi), 2.*M_PI*x[k++]);
		eij = std::polar(xi, 2.*M_PI*x[k++]);
		coli = M.col(i);
		M.col(i) = M.col(i)*eii - M.col(j)*conj(eij);
		M.col(j) =     coli*eij + M.col(j)*conj(eii);
	}
	return;
}

matrix_type linopt::hurwitz_parametrization(const point &x)
{
	matrix_type M;
	hurwitz_parametrization(M, x);
	return M;
}

void linopt::exp_hermite_parametrization(matrix_type &M, const point &x)
{
	const int N = std::sqrt(x.size());
	assert(N*N == static_cast<int>(x.size()));
	matrix_type H(N, N);
	int k = N;
	for(int i = 0; i < N; i++)
	{
		H(i, i) = x[i];
		for(int j = i + 1; j < N; j++)
		{
			H(i, j) = complex_type(x[k],  x[k+1]);
			H(j, i) = complex_type(x[k], -x[k+1]);
			k += 2;
		}
	}
	H *= complex_type(0, 1);
	M = H.exp();
	return;
}

matrix_type linopt::exp_hermite_parametrization(const point &x)
{
	matrix_type M;
	exp_hermite_parametrization(M, x);
	return M;
}

bool linopt::is_column_unitary(const matrix_type &M, real_type eps)
{
	return (M.adjoint() * M).isIdentity(eps);
}

bool linopt::is_row_unitary(const matrix_type &M, real_type eps)
{
	return (M * M.adjoint()).isIdentity(eps);
}

bool linopt::is_unitary(const matrix_type &M, real_type eps)
{
	if(M.cols() != M.rows())
		return false;
	else
		return is_column_unitary(M, eps) && is_row_unitary(M, eps);
}
