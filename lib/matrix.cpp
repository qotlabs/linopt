// SPDX-License-Identifier: LGPL-3.0-or-later
// SPDX-FileCopyrightText: Copyright (c) 2018-2025, Quantum Optical Technologies Laboratories
// SPDX-FileContributor: Struchalin Gleb <struchalin.gleb@physics.msu.ru>
// SPDX-FileContributor: Dyakonov Ivan <iv.dyakonov@physics.msu.ru>

#include "matrix.h"
#include "misc.h"
#include "exceptions.h"
#include <limits.h>
#include <cmath>
#include <complex>
#include <unsupported/Eigen/MatrixFunctions>

using namespace linopt;

Complex linopt::permanent(const Matrix &M)
{
	if(M.cols() != M.rows())
		throw WrongSize(ERROR_MSG("Matrix should be square. Current size is (" +
				std::to_string(M.rows()) + "x" + std::to_string(M.cols()) + ")."));
	if(M.cols() == 0)
		throw WrongSize(ERROR_MSG("Empty matrix is not allowed."));
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
	using Index = unsigned long;
	if(static_cast<Index>(M.cols()) > CHAR_BIT*sizeof(Index))
		throw WrongSize(ERROR_MSG("Matrix size (which is " +
				std::to_string(M.rows()) + ") is too large (maximum is " +
				std::to_string(CHAR_BIT*sizeof(Index)) + ")."));
	Complex perm = 0.;
	Column sum = M.rowwise().sum();
	Real mult;
	Index n = 0;
	Index nmax = static_cast<Index>(1) << (M.cols() - 1);
	while(n < nmax)
	{
		Complex sp = sum.prod();
		(n & 1) ? perm -= sp : perm += sp;
		Index i = (n ^ (n + 1)) + 1;
		n++;
		mult = (n & i) ? 2 : -2;
		i = ctz(i) - 1;
		sum += mult * M.col(i);
	}
	perm /= nmax;
	return perm;
}

static inline Real pyramid(Real x, Real a)
{
	x = std::fabs(x);
	x = ((unsigned)(x/a) & 1) ? -x : x;
	return mod(x, a);
}

void linopt::hurwitzParametrization(Matrix &M, const Point &x)
{
	const int N = isqrt(x.size());
	if(N*N != static_cast<int>(x.size()))
		throw WrongSize(ERROR_MSG("Size of point (which is " +
				std::to_string(x.size()) + ") should be a square of an integer."));
	if(M.cols() != N || M.rows() != N)
		M.resize(N, N);
	int i, j, k = 0;
	Complex eii, eij;
	Column coli(N);
	Real xi;
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

Matrix linopt::hurwitzParametrization(const Point &x)
{
	Matrix M;
	hurwitzParametrization(M, x);
	return M;
}

void linopt::expHermiteParametrization(Matrix &M, const Point &x)
{
	const int N = isqrt(x.size());
	if(N*N != static_cast<int>(x.size()))
		throw WrongSize(ERROR_MSG("Size of point (which is " +
				std::to_string(x.size()) + ") should be a square of an integer."));
	Matrix H(N, N);
	int k = N;
	for(int i = 0; i < N; i++)
	{
		H(i, i) = x[i];
		for(int j = i + 1; j < N; j++)
		{
			H(i, j) = Complex(x[k],  x[k+1]);
			H(j, i) = Complex(x[k], -x[k+1]);
			k += 2;
		}
	}
	H *= Complex(0, 1);
	M = H.exp();
	return;
}

Matrix linopt::expHermiteParametrization(const Point &x)
{
	Matrix M;
	expHermiteParametrization(M, x);
	return M;
}

Real linopt::matrixFidelity(const Matrix &A, const Matrix &B)
{
	if(A.rows() != B.rows() || A.cols() != B.cols())
		throw WrongSize(ERROR_MSG("Matrix A size (which is " +
			std::to_string(A.rows()) + "x" + std::to_string(A.cols()) +
			") should be equal to B size (which is " +
			std::to_string(B.rows()) + "x" + std::to_string(B.cols()) + ")."));
	Real F = std::norm((A.conjugate().cwiseProduct(B)).sum());
	F /= A.squaredNorm() * B.squaredNorm();
	return F;
}

bool linopt::isColumnUnitary(const Matrix &M, Real eps)
{
	return (M.adjoint() * M).isIdentity(eps);
}

bool linopt::isRowUnitary(const Matrix &M, Real eps)
{
	return (M * M.adjoint()).isIdentity(eps);
}

bool linopt::isUnitary(const Matrix &M, Real eps)
{
	if(M.cols() != M.rows())
		return false;
	else
		return isColumnUnitary(M, eps) && isRowUnitary(M, eps);
}
