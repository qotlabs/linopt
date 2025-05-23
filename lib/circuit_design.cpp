// SPDX-License-Identifier: LGPL-3.0-or-later
// SPDX-FileCopyrightText: Copyright (c) 2018-2025, Quantum Optical Technologies Laboratories
// SPDX-FileContributor: Struchalin Gleb <struchalin.gleb@physics.msu.ru>
// SPDX-FileContributor: Dyakonov Ivan <iv.dyakonov@physics.msu.ru>
// SPDX-FileContributor: Fldjyan Suren <fldzhian.sa17@physics.msu.ru>

#include "circuit_design.h"
#include "exceptions.h"
#include <memory>

using namespace linopt;

static void checkmate(int mm[], const int N)
{
	// mm is assumed to contain N(N-1)/2 elements
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

void linopt::clementsDesign(Matrix &M, const Point &x, const Point &y)
{
	// Gives unitary matrix assured
	if(x.size() != y.size())
		throw WrongSize(ERROR_MSG("Size of x (which is " + std::to_string(x.size()) +
			") should be equal to size of y (which is " + std::to_string(y.size()) +
			")."));
	const int NN = static_cast<int>(x.size() / 2.);
	const int N = static_cast<int>(std::sqrt(2.*NN + 0.25) + 0.5);
	if(N*(N - 1) / 2 != NN)
		throw WrongSize(ERROR_MSG("Wrong size of points (current is " +
			std::to_string(x.size()) + ")."));
	if(M.cols() != N || M.rows() != N || !isUnitary(M))
	{
		M.resize(N, N);
		M.setIdentity();
	}
	std::unique_ptr<int[]> mm(new int[NN]);
	checkmate(mm.get(), N);
	Column colm(N);
	const Complex i1(0., 1.);
	Complex expphi, exp2theta;
	Real cosp, cosm, sinp, sinm;
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
	return;
}

void linopt::clementsDesign(Matrix &M, const Point &x)
{
	// Gives unitary matrix assured
	const int NN = static_cast<int>(x.size() / 2.);
	const int N = static_cast<int>(std::sqrt(2.*NN + 0.25) + 0.5);
	if(N*(N - 1) / 2 != NN)
		throw WrongSize(ERROR_MSG("Wrong size of point (current is " +
			std::to_string(x.size()) + ")."));
	if (M.cols() != N || M.rows() != N || !isUnitary(M))
	{
		M.resize(N, N);
		M.setIdentity();
	}
	std::unique_ptr<int[]> mm(new int[NN]);
	checkmate(mm.get(), N);
	Column colm(N);
	const Complex i1(0., 1.);
	Complex expphi, exp2theta;
	for(int i = 0; i < NN; i++)
	{
		int m = mm[NN - 1 - i] - 1;
		expphi = std::polar(1., x[2 * i]);
		exp2theta = std::polar(1., 2.*x[2 * i + 1]);
		colm = M.col(m);
		M.col(m)   = M.col(m)*((exp2theta - 1.)*expphi*0.5) + M.col(m + 1)*((exp2theta + 1.)*i1*expphi*0.5);
		M.col(m+1) = colm    *((exp2theta + 1.)*i1*0.5)     - M.col(m + 1)*((exp2theta - 1.)*0.5);
	}
	return;
}

Matrix linopt::clementsDesign(const Point &x, const Point &y)
{
	Matrix M;
	clementsDesign(M, x, y);
	return M;
}

Matrix linopt::clementsDesign(const Point &x)
{
	Matrix M;
	clementsDesign(M, x);
	return M;
}

void linopt::getClementsDesign(Matrix &M, Point &x, Real eps)
{
	// Leaves diagonal unitary matrix & gives point in which even elements are
	// phi, odd elements are theta. Phase shift coefficients are given in
	// reverse column-wise enumeration order
	const int N = static_cast<int>(M.cols());
	if(eps > 0)
		if (!isUnitary(M, eps))
			throw NotUnitary(ERROR_MSG("Given matrix is not unitary"));
	const int NN = N * (N - 1) / 2;
	if(static_cast<int>(x.size()) != 2*NN)
		x.resize(2*NN);
	std::unique_ptr<int[]> mm(new int[NN]);
	std::unique_ptr<int[]> pattern(new int[NN]);
	checkmate(pattern.get(), N);
	int f = NN - 1, s = 0, i, n, m;
	const Complex i1(0., 1.);
	Complex expphi, exp2theta;
	Real phi, theta, phim, phin;
	Column tempcol(N);
	Row temprow(N);
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
				x[2*f] = phi;
				x[2*f + 1] = theta;
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
				x[2 * s] = phi;
				x[2 * s + 1] = theta;
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
		M(m, m) = -polar(1., phin - x[2 * i]);
		if (phim - phin >= 0)
			x[2 * i] = phim - phin - M_PI;
		else
			x[2 * i] = phim - phin + M_PI;
		x[2 * i + 1] = -x[2 * i + 1];
	}
	// Obtaining reverse column-wise enumeration order
	for(i = NN - 1; i >= 0; i--)
	{
		if(mm[i] != pattern[NN - 1 - i])
		{
			for(n = i - 1; n >= 0; n--)
			{
				if(mm[n] == pattern[NN - 1 - i])
				{
					phim = x[2 * i];
					x[2 * i] = x[2 * n];
					x[2 * n] = phim;
					phim = x[2 * i + 1];
					x[2 * i + 1] = x[2 * n + 1];
					x[2 * n + 1] = phim;
					f = mm[i];
					mm[i] = mm[n];
					mm[n] = f;
					break;
				}
				else if(mm[n] == mm[i])
				{
					phim = x[2 * i];
					x[2 * i] = x[2 * n];
					x[2 * n] = phim;
					phim = x[2 * i + 1];
					x[2 * i + 1] = x[2 * n + 1];
					x[2 * n + 1] = phim;
				}
			}
		}
	}
}

Point linopt::getClementsDesign(Matrix &M, Real eps)
{
	Point x;
	getClementsDesign(M, x, eps);
	return x;
}
