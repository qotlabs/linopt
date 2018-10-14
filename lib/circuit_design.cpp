/* Copyright © 2018, Quantum Optical Technologies Laboratories
 * <https://www.qotlabs.org/en/>
 * Contributed by: Struchalin Gleb <struchalin.gleb@physics.msu.ru>
 *                 Dyakonov Ivan <iv.dyakonov@physics.msu.ru>
 *				   Fldjyan Suren
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

/** @defgroup design Circuit design
 * @brief Circuit designs and decompositions into elementary blocks
 */

#include "circuit_design.h"
#include "exceptions.h"

using namespace linopt;

static void checkmate(int mm[], const int N)
{
	// mm is assumed to contain N(N-1)/2 elementss
	int k = 1;
	for(int i = 0; i < N*(N - 1) / 2; i++)
	{
		mm[i] = k;
		k += 2;
		if(k > N - 1)
			k = mm[i] % 2 + 1;
	}
	return;
}

/** @ingroup design
 * @brief
 *
 * @param M[out] --
 * @param x[in] --
 * @param y[in] --
 *
 * @throw
 * If `x` size is not equal to @f$ N(N-1)/2 @f$ for some integer @f$ N @f$ or if
 * `x` and `y` sizes do not match, then `wrong_size` is thrown.
 *
 * @see
 * W. R. Clements, _et al._ "Optimal design for universal multiport
 * interferometers." Optica __3__, pp. 1460-1465 (2016),
 * https://doi.org/10.1364/OPTICA.3.001460
 */
void linopt::clements_design(matrix_type &M, const point &x, const point &y)
{
	// Gives unitary matrix assured
	if(x.size() != y.size())
		throw wrong_size(ERROR_MSG("Size of x (which is " + std::to_string(x.size()) +
			") should be equal to size of y (which is " + std::to_string(y.size()) +
			")."));
	const int NN = static_cast<int>(x.size() / 2.);
	const int N = static_cast<int>(std::sqrt(2.*NN + 0.25) + 0.5);
	if(N*(N - 1) / 2 != NN)
		throw wrong_size(ERROR_MSG("Wrong size of points (current is " +
			std::to_string(x.size()) + ")."));
	if(M.cols() != N || M.rows() != N)
	{
		M.resize(N, N);
		M.setIdentity();
	}
	int* mm = new int[NN];
	checkmate(mm, N);
	matrix_type colm(N, 1);
	const complex_type i1(0., 1.);
	complex_type expphi, exp2theta;
	real_type cosp, cosm, sinp, sinm;
	using std::cos;
	using std::sin;
	using std::polar;
	for(int i = 0; i < NN; i++)
	{
		int m = mm[NN - 1 - i] - 1;
		cosp = cos(y[2 * i] + y[2 * i + 1]);
		cosm = cos(y[2 * i] - y[2 * i + 1]);
		sinp = sin(y[2 * i] + y[2 * i + 1]);
		sinm = sin(y[2 * i] - y[2 * i + 1]);
		expphi = polar(1., x[2 * i]);
		exp2theta = polar(1., 2.*x[2 * i + 1]);
		colm = M.col(m);
		M.col(m)   = M.col(m)  *((exp2theta*(cosm - sinp) - (cosm + sinp))*expphi*0.5) +
					 M.col(m+1)*((exp2theta*(cosp - sinm) + (cosp + sinm))*i1*expphi*0.5);
		M.col(m+1) = colm      *((exp2theta*(cosp + sinm) + (cosp - sinm))*i1*0.5) -
					 M.col(m+1)*((exp2theta*(cosm + sinp) - (cosm - sinp))*0.5);
	}
	delete[] mm;
	return;
}

/** @ingroup design
 * @brief
 *
 * @overload
 */
void linopt::clements_design(matrix_type &M, const point &x)
{
	// Gives unitary matrix assured
	const int NN = static_cast<int>(x.size() / 2.);
	const int N = static_cast<int>(std::sqrt(2.*NN + 0.25) + 0.5);
	if(N*(N - 1) / 2 != NN)
		throw wrong_size(ERROR_MSG("Wrong size of point (current is " +
			std::to_string(x.size()) + ")."));
	if(M.cols() != N || M.rows() != N)
	{
		M.resize(N, N);
		M.setIdentity();
	}
	int *mm = new int[NN];
	checkmate(mm, N);
	matrix_type colm(N, 1);
	const complex_type i1(0., 1.);
	complex_type expphi, exp2theta;
	for(int i = 0; i < NN; i++)
	{
		int m = mm[NN - 1 - i] - 1;
		expphi = std::polar(1., x[2 * i]);
		exp2theta = std::polar(1., 2.*x[2 * i + 1]);
		colm = M.col(m);
		M.col(m)   = M.col(m)*((exp2theta - 1.)*expphi*0.5) + M.col(m + 1)*((exp2theta + 1.)*i1*expphi*0.5);
		M.col(m+1) = colm    *((exp2theta + 1.)*i1*0.5)     - M.col(m + 1)*((exp2theta - 1.)*0.5);
	}
	delete[] mm;
	return;
}

matrix_type linopt::clements_design(const point &x, const point &y)
{
	matrix_type M;
	clements_design(M, x, y);
	return M;
}

matrix_type linopt::clements_design(const point &x)
{
	matrix_type M;
	clements_design(M, x);
	return M;
}

point linopt::get_clements_design(matrix_type &M)
{
	// Leaves diagonal unitary matrix & gives point in which even elements are
	// phi, odd elements are theta. Angle coefficients are given in
	// multiplicative order
	const int N = static_cast<int>(M.cols());
	if (M.cols() != N || M.rows() != N)
		M.resize(N, N);
	const int NN = N * (N - 1) / 2;
	point param(2 * NN);
	int* mm = new int[NN];
	int* pattern = new int[NN];
	checkmate(pattern, N);
	int f = NN - 1, s = 0, i, n, m;
	const complex_type i1(0., 1.);
	complex_type expphi, exp2theta;
	real_type phi, theta, phim, phin;
	matrix_type tempcol(N, 1), temprow(1, N);
	using std::arg;
	using std::atan;
	using std::abs;
	using std::polar;
	// Matrix diagonalization
	for(int k = 0; k < N - 1; k++)
	{
		if(k % 2 == 0)
		{
			for(m = k; m >= 0; m--)
			{
				i = N - 1 - k + m;
				mm[f] = m + 1;
				phi = -arg(-M(i, m + 1) / M(i, m));
				theta = atan(abs(M(i, m + 1) / M(i, m)));
				tempcol = M.col(m);
				expphi = polar(1., -phi);
				exp2theta = polar(1., -2. * theta);
				M.col(m)   = M.col(m)  *((exp2theta - 1.)*expphi*0.5) -
							 M.col(m+1)*((exp2theta + 1.)*i1*0.5);
				M.col(m+1) = tempcol   *((exp2theta + 1.)*i1*expphi*(-0.5)) -
							 M.col(m+1)*((exp2theta - 1.)*0.5);
				param[2*f] = phi;
				param[2*f + 1] = theta;
				f--;
			}
		}
		else
		{
			for(int j = 0; j <= k; j++)
			{
				n = N - 1 - k + j;
				mm[s] = n;
				phi = -arg(M(n - 1, j) / M(n, j));
				theta = atan(abs(M(n - 1, j) / M(n, j)));
				temprow = M.row(n - 1);
				expphi = polar(1., phi);
				exp2theta = polar(1., 2. * theta);
				M.row(n-1) = M.row(n-1)*((exp2theta - 1.)*expphi*0.5) +
							 M.row(n)  *((exp2theta + 1.)*i1*0.5);
				M.row(n)   = temprow   *((exp2theta + 1.)*expphi*i1*0.5) -
							 M.row(n)  *((exp2theta - 1.)*0.5);
				param[2 * s] = phi;
				param[2 * s + 1] = theta;
				s++;
			}
		}
	}
	// Moving ancillary matrices
	s--;
	for(i = s; i >= 0; i--)
	{
		m = mm[i] - 1;
		phim = arg(M(m, m));
		phin = arg(M(m + 1, m + 1));
		M(m, m) = -polar(1., phin - param[2 * i]);
		if (phim - phin >= 0)
			param[2 * i] = phim - phin - M_PI;
		else
			param[2 * i] = phim - phin + M_PI;
		param[2 * i + 1] = -param[2 * i + 1];
	}
	// Obtaining checkmate-friendly order
	for(i = NN - 1; i >= 0; i--)
	{
		if(mm[i] != pattern[NN - 1 - i])
		{
			for(n = i - 1; n >= 0; n--)
			{
				if(mm[n] == pattern[NN - 1 - i])
				{
					phim = param[2 * i];
					param[2 * i] = param[2 * n];
					param[2 * n] = phim;
					phim = param[2 * i + 1];
					param[2 * i + 1] = param[2 * n + 1];
					param[2 * n + 1] = phim;
					f = mm[i];
					mm[i] = mm[n];
					mm[n] = f;
					break;
				}
				else if(mm[n] == mm[i])
				{
					phim = param[2 * i];
					param[2 * i] = param[2 * n];
					param[2 * n] = phim;
					phim = param[2 * i + 1];
					param[2 * i + 1] = param[2 * n + 1];
					param[2 * n + 1] = phim;
				}
			}
		}
	}
	delete[] mm;
	delete[] pattern;
	return param;
}
