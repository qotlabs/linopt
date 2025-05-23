// SPDX-License-Identifier: LGPL-3.0-or-later
// SPDX-FileCopyrightText: Copyright (c) 2018-2025, Quantum Optical Technologies Laboratories
// SPDX-FileContributor: Struchalin Gleb <struchalin.gleb@physics.msu.ru>
// SPDX-FileContributor: Dyakonov Ivan <iv.dyakonov@physics.msu.ru>

#pragma once

#include "types.h"
#include <Eigen/Dense>
#include <vector>

/** @defgroup matrix Matrix
 * @brief Various matrix-related functions.
 */

namespace linopt
{

using Matrix = Eigen::Matrix<Complex, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>;
using Column = Eigen::Matrix<Complex, Eigen::Dynamic, 1>;
using Row = Eigen::Matrix<Complex, 1, Eigen::Dynamic>;
using Point = std::vector<Real>;

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
 * `WrongSize` is thrown.
 *
 * @see
 * K. Życzkowski and M. Kus "Random unitary matrices." J. Phys. A: Math. Gen.
 * __27__, pp. 4235-4245 (1994), https://doi.org/10.1088/0305-4470/27/12/028
 */
void hurwitzParametrization(Matrix &M, const Point &x);

/** @ingroup matrix
 * @brief Returns a matrix parametrized according to the Hurwitz
 * parametrization.
 *
 * @param[in] x -- an array of parameters.
 * @return Unitary matrix parametrized by `x`.
 *
 * @overload
 */
Matrix hurwitzParametrization(const Point &x);

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
 * `WrongSize` is thrown.
 */
void expHermiteParametrization(Matrix &M, const Point &x);

/** @ingroup matrix
 * @brief Returns a matrix parametrized according to the exponential-Hermitian
 * parametrization.
 *
 * @param[in] x -- an array of parameters.
 * @return Unitary matrix parametrized by `x`.
 *
 * @overload
 */
Matrix expHermiteParametrization(const Point &x);

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
 * @throw If matrix sizes do not match `WrongSize` is thrown.
 */
Real matrixFidelity(const Matrix &A, const Matrix &B);

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
 * @see isRowUnitary(), isUnitary()
 */
bool isColumnUnitary(const Matrix &M, Real eps = defaultEpsilon);

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
 * @see isColumnUnitary(), isUnitary()
 */
bool isRowUnitary(const Matrix &M, Real eps = defaultEpsilon);

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
 * @see isColumnUnitary(), isRowUnitary()
 */
bool isUnitary(const Matrix &M, Real eps = defaultEpsilon);

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
Complex permanent(const Matrix &M);

} // namespace linopt
