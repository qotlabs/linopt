/* Copyright © 2018, Quantum Optical Technologies Laboratories
 * <https://www.qotlabs.org/en/>
 * Contributed by: Struchalin Gleb <struchalin.gleb@physics.msu.ru>
 *                 Dyakonov Ivan <iv.dyakonov@physics.msu.ru>
 *
 * This file is part of Linopt.
 *
 * Linopt is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Linopt is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with Linopt. If not, see <https://www.gnu.org/licenses/>.
 */

/** @defgroup matrix Matrix
 * @brief Various matrix-related functions.
 */

#include <limits.h>
#include <cmath>
#include <complex>
#include <unsupported/Eigen/MatrixFunctions>
#include "misc.h"
#include "matrix.h"
#include "exceptions.h"

using namespace linopt;

/** @ingroup matrix
 * @brief Calculates permanent of a square matrix.
 *
 * @param[in] M -- matrix to calculate its permanent.
 * @return Computed permanent.
 *
 * Calculates permanent of a square @f$ N \times N @f$ matrix `M`.
 * Internally the function uses Glynn formula with Gray code summation
 * technique. Complexity is @f$ O(N \times 2^N) @f$.
 *
 * @todo Parallel version.
 *
 * @see
 * David G.Glynn "The permanent of a square matrix." Eur. J. Combin. __31__,
 * pp. 1887-1891 (2010), https://doi.org/10.1016/j.ejc.2010.01.010
 */
complex_type linopt::permanent(const matrix_type &M)
{
	if(M.cols() != M.rows())
		throw wrong_size(ERROR_MSG("Matrix should be square. Current size is (" +
				std::to_string(M.rows()) + "x" + std::to_string(M.cols()) + ")."));
	if(M.cols() == 0)
		throw wrong_size(ERROR_MSG("Empty matrix is not allowed."));
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
	using index_type = unsigned long;
	if(static_cast<index_type>(M.cols()) > CHAR_BIT*sizeof(index_type))
		throw wrong_size(ERROR_MSG("Matrix size (which is " +
				std::to_string(M.rows()) + ") is too large (maximum is " +
				std::to_string(CHAR_BIT*sizeof(index_type)) + ")."));
	complex_type perm = 0.;
	column_type sum = M.rowwise().sum();
	real_type mult;
	index_type n = 0;
	index_type nmax = static_cast<index_type>(1) << (M.cols() - 1);
	while(n < nmax)
	{
		complex_type sp = sum.prod();
		(n & 1) ? perm -= sp : perm += sp;
		index_type i = (n ^ (n + 1)) + 1;
		n++;
		mult = (n & i) ? 2 : -2;
		i = ctz(i) - 1;
		sum += mult * M.col(i);
	}
	perm /= nmax;
	return perm;
}

static inline real_type pyramid(real_type x, real_type a)
{
	x = std::fabs(x);
	x = ((unsigned)(x/a) & 1) ? -x : x;
	return mod(x, a);
}

/** @ingroup matrix
 * @brief Fills a unitary matrix according to the Hurwitz parametrization.
 *
 * @param[out] M -- an @f$ N \times N @f$ matrix to store the result. If it has
 * an improper size it will be resized.
 * @param[in] x -- an array of @f$ N^2 @f$ parameters.
 *
 * Hurwitz parametrization is capable to parametrize any unitary matrix of size
 * @f$ N \times N @f$ with @f$ N^2 @f$ parameters `x`. The essential feature of
 * this parametrization is that a uniform ensemble of random points inside
 * a unit hypercube @f$ 0 \le x_i < 1, i = 1, \dots, N^2 @f$ is mapped onto
 * a Haar-uniform ensemble of random unitary matrices. Therefore this
 * parametrization can be used, for example, to sample random matrices from
 * a circular unitary ensemble (CUE). Or it can be used when optimizing over
 * the space of all unitary matrices. In some sense it maps the space of unitary
 * matrices onto a unit hypercube with least possible "distorsion".
 *
 * Parameters `x` may lie outside the unit hypercube. If so, they are
 * effectively moved back to the hypercube by means of a peridic function.
 * The parametrization is continuous on boundaries, but not smooth.
 *
 * @throw
 * If size of `x` is not a square of some integer number @f$ N @f$, then
 * `wrong_size` is thrown.
 *
 * @see
 * K. Życzkowski and M. Kus "Random unitary matrices." J. Phys. A: Math. Gen.
 * __27__, pp. 4235-4245 (1994), https://doi.org/10.1088/0305-4470/27/12/028
 */
void linopt::hurwitz_parametrization(matrix_type &M, const point &x)
{
	const int N = isqrt(x.size());
	if(N*N != static_cast<int>(x.size()))
		throw wrong_size(ERROR_MSG("Size of point (which is " +
				std::to_string(x.size()) + ") should be a square of an integer."));
	if(M.cols() != N || M.rows() != N)
		M.resize(N, N);
	int i, j, k = 0;
	complex_type eii, eij;
	column_type coli(N);
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

/** @ingroup matrix
 * @brief Returns a matrix parametrized according to the Hurwitz
 * parametrization.
 *
 * @param[in] x -- an array of parameters.
 * @return Unitary matrix parametrized by `x`.
 *
 * @overload
 */
matrix_type linopt::hurwitz_parametrization(const point &x)
{
	matrix_type M;
	hurwitz_parametrization(M, x);
	return M;
}

/** @ingroup matrix
 * @brief Fills a unitary matrix according to the exponential-Hermitian
 * parametrization.
 *
 * @param[out] M -- an @f$ N \times N @f$  matrix to store the result. If it has
 * improper size it will be resized.
 * @param[in] x -- an array of @f$ N^2 @f$ parameters.
 *
 * The exponential-Hermitian parametrization of a unitary matrix is based on
 * the relation
 * @f[ M = \exp(i H(x)). @f]
 * If @f$ H @f$ is Hermitian then @f$ M @f$ is unitary. Since @f$ H @f$ is
 * Hermitian it can be parametrized by its matrix elements. First elements of
 * `x` contain the real-valued diagonal of @f$ H @f$, then the upper triangular
 * part is specified. For example, @f$ 4 \times 4 @f$ Hermitian matrix is
 * parametrized in the following way
 * @f[
 *	H(x) =
 *	\begin{pmatrix}
 *	    x[0]    &  x[4]+ix[5]  &  x[6]+ix[7]  &  x[8]+ix[9]  \\
 *	x[4]-ix[5]  &      x[1]    & x[10]+ix[11] & x[12]+ix[13] \\
 *	x[6]-ix[7]  & x[10]-ix[11] &     x[2]     & x[14]+ix[15] \\
 *	x[8]-ix[9]  & x[12]-ix[13] & x[14]-ix[15] &      x[3]
 *	\end{pmatrix}.
 * \f]
 *
 * @throw
 * If size of `x` is not a square of some integer number @f$ N @f$, then
 * `wrong_size` is thrown.
 */
void linopt::exp_hermite_parametrization(matrix_type &M, const point &x)
{
	const int N = isqrt(x.size());
	if(N*N != static_cast<int>(x.size()))
		throw wrong_size(ERROR_MSG("Size of point (which is " +
				std::to_string(x.size()) + ") should be a square of an integer."));
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

/** @ingroup matrix
 * @brief Returns a matrix parametrized according to the exponential-Hermitian
 * parametrization.
 *
 * @param[in] x -- an array of parameters.
 * @return Unitary matrix parametrized by `x`.
 *
 * @overload
 */
matrix_type linopt::exp_hermite_parametrization(const point &x)
{
	matrix_type M;
	exp_hermite_parametrization(M, x);
	return M;
}

/** @ingroup matrix
 * @brief Calculates normalized fidelity between two matrices.
 *
 * @param[in] A -- first matrix.
 * @param[in] B -- second matrix.
 * @return Calculated fidelity.
 *
 * Normalized fidelity @f$ F(A, B) @f$ between two matrices @f$ A @f$ and
 * @f$ B @f$ is defined as follows:
 * @f[ \newcommand{\Tr}{\mathrm{Tr}}
 *     F(A, B) = \frac{|\Tr(A^\dagger B)|^2}{\Tr(A^\dagger A) \Tr(B^\dagger B)}.
 * @f]
 *
 * @throw If matrix sizes do not match `wrong_size` is thrown.
 */
real_type linopt::matrix_fidelity(const matrix_type &A, const matrix_type &B)
{
	if(A.rows() != B.rows() || A.cols() != B.cols())
		throw wrong_size(ERROR_MSG("Matrix A size (which is " +
			std::to_string(A.rows()) + "x" + std::to_string(A.cols()) +
			") should be equal to B size (which is " +
			std::to_string(B.rows()) + "x" + std::to_string(B.cols()) + ")."));
	real_type F = std::norm((A.conjugate().cwiseProduct(B)).sum());
	F /= A.squaredNorm() * B.squaredNorm();
	return F;
}

/** @ingroup matrix
 * @brief Tests whether a matrix is column-unitary.
 *
 * @param[in] M -- matrix to test.
 * @param[in] eps -- accuracy within test is performed.
 * @return `true` if the matrix `M` is a column-unitary one with the given
 * accuracy, `false` otherwise.
 *
 * A matrix @f$ M @f$ is called column unitary if @f$ M^\dagger M = 1 @f$.
 *
 * @see is_row_unitary(), is_unitary()
 */
bool linopt::is_column_unitary(const matrix_type &M, real_type eps)
{
	return (M.adjoint() * M).isIdentity(eps);
}

/** @ingroup matrix
 * @brief Tests whether a matrix is row-unitary.
 *
 * @param[in] M -- matrix to test.
 * @param[in] eps -- accuracy within test is performed.
 * @return `true` if the matrix `M` is a row-unitary one with the given
 * accuracy, `false` otherwise.
 *
 * A matrix @f$ M @f$ is called row unitary if @f$ M M^\dagger = 1 @f$.
 *
 * @see is_column_unitary(), is_unitary()
 */
bool linopt::is_row_unitary(const matrix_type &M, real_type eps)
{
	return (M * M.adjoint()).isIdentity(eps);
}

/** @ingroup matrix
 * @brief Tests whether a matrix is unitary.
 *
 * @param[in] M -- matrix to test.
 * @param[in] eps -- accuracy within test is performed.
 * @return `true` if the matrix `M` is a unitary one with the given
 * accuracy, `false` otherwise.
 *
 * A square matrix @f$ M @f$ is called unitary if it is both column and row
 * unitary.
 *
 * @see is_column_unitary(), is_row_unitary()
 */
bool linopt::is_unitary(const matrix_type &M, real_type eps)
{
	if(M.cols() != M.rows())
		return false;
	else
		return is_column_unitary(M, eps) && is_row_unitary(M, eps);
}
