// SPDX-License-Identifier: LGPL-3.0-or-later
// SPDX-FileCopyrightText: Copyright (c) 2018-2025, Quantum Optical Technologies Laboratories
// SPDX-FileContributor: Struchalin Gleb <struchalin.gleb@physics.msu.ru>
// SPDX-FileContributor: Dyakonov Ivan <iv.dyakonov@physics.msu.ru>
// SPDX-FileContributor: Fldjyan Suren <fldzhian.sa17@physics.msu.ru>

#pragma once

#include "matrix.h"

/** @defgroup design Circuit design
 * @brief Circuit designs and decompositions into elementary blocks.
 */

namespace linopt
{

/** @ingroup design
 * @brief Fills a unitary matrix according to the Clements design.
 *
 * @param[in,out] M -- on input a diagonal unitary matrix of size
 * @f$ N \times N @f$ is provided. If the matrix is not unitary or have improper
 * size (calculated from `x` size), then it will be resized accordingly and
 * set to the identity matrix. On output this matrix is used to store the final
 * result.
 * @param[in] x -- an array of @f$ N(N-1) @f$ phase-shift parameters, such
 * that even elements are @f$ \phi @f$ -- phase shifts _before_ the beam
 * splitters, and odd ones are @f$ \theta @f$ -- halves of phase shifts _between_ the beam
 * splitters. All pairs of parameters go in reverse column-wise enumeration
 * order.
 * @param[in] y -- an array of beam splitters angle defects, such that even
 * elements are defects of the first splitters, and odd ones are defects of the
 * second splitters. All pairs of parameters go in reverse column-wise
 * enumeration order.
 *
 * Fills a unitary matrix parametrized according to the Clements parametrization
 * with basic matrices of the form
 * @f[
 *		i e^{i\theta} \begin{pmatrix}
 *		e^{i\phi}\sin{\theta} &  \cos{\theta}  \\
 *		e^{i\phi}\cos{\theta} & -\sin{\theta}
 *  \end{pmatrix}
 * @f]
 * with beam splitter angle defects.
 *
 * @throw
 * If `x` size is not equal to @f$ N(N-1) @f$ for some integer @f$ N @f$ or if
 * `x` and `y` sizes do not match, then `WrongSize` is thrown.
 *
 * @see getClementsDesign()
 *
 * @see
 * W. R. Clements, _et al._ "Optimal design for universal multiport
 * interferometers." Optica __3__, pp. 1460-1465 (2016),
 * https://doi.org/10.1364/OPTICA.3.001460
 */
void clementsDesign(Matrix &M, const Point &x, const Point &y);

/** @ingroup design
 * @brief Equivalent to `clementsDesign(M, x, y)` with all @f$ y_i = 0 @f$.
 *
 * This function is equivalent to the `clementsDesign(M, x, y)` with all
 * elements of `y` set to zero. However, calculations are a bit faster than
 * direct call of `clementsDesign(M, x, y)`.
 */
void clementsDesign(Matrix &M, const Point &x);

/** @ingroup design
 * @brief This function is equivalent to `clementsDesign(M, x, y)` with `M`
 * being the identity matrix.
 *
 * @return Calculated unitary matrix `M`.
 */
Matrix clementsDesign(const Point &x, const Point &y);

/** @ingroup design
 * @brief This function is equivalent to `clementsDesign(M, x)` with `M` being
 * the identity matrix.
 *
 * @return Calculated unitary matrix `M`.
 */
Matrix clementsDesign(const Point &x);

/** @ingroup design
 * @brief Calculates phase-shift coefficients for a unitary matrix `M` according
 * to the Clements design.
 *
 * @param[in,out] M -- on input a unitary @f$ N \times N @f$ matrix is
 * specified. On output a diagonal unitary matrix is produced. The original
 * matrix is destroyed during calculations.
 * @param[out] x -- the array of @f$ N(N-1) @f$ phase-shift parameters, such
 * that even elements are @f$ \phi @f$ -- phase shifts _before_ the beam
 * splitters, and odd ones are @f$ \theta @f$ -- halves of phase shifts
 * _between_ the beam splitters. All pairs of parameters go in reverse
 * column-wise enumeration order. If `x` has improper size then it will be
 * resized.
 * @param[in] eps -- precision for unitarity test of the input matrix `M`. If
 * `eps` is negative then no tests are performed.
 *
 * Returns phase-shift coefficients and leaves a diagonal unitary matrix
 * according to the Clements design with basic matrices of the form
 * @f[
 *		i e^{i\theta} \begin{pmatrix}
 *		e^{i\phi}\sin{\theta} &  \cos{\theta}  \\
 *		e^{i\phi}\cos{\theta} & -\sin{\theta}
 *	\end{pmatrix}.
 * @f]
 * As a result of multiplication
 * @f[
 *
 *	\begin{pmatrix}
 *		1/\sqrt{2} &  i/\sqrt{2}  \\
 *		i/\sqrt{2} &  1/\sqrt{2}
 *	\end{pmatrix}
 *	\begin{pmatrix}
 *		e^{i2\theta} &  0  \\
 *		0 &  1
 *	\end{pmatrix}
 *	\begin{pmatrix}
 *		1/\sqrt{2} &  i/\sqrt{2}  \\
 *		i/\sqrt{2} &  1/\sqrt{2}
 * 	\end{pmatrix}
 *	\begin{pmatrix}
 *		e^{i\phi} &  0  \\
 *		0 &  1
 *		\end{pmatrix}
 *
 * @f]
 * This function is effectively inverse of the `clementsDesign(M, x)`.
 *
 * @throw
 * If `M` is not unitary within given precision `eps`, then `NotUnitary` is
 * thrown.
 *
 * @see clementsDesign()
 *
 * @see
 * W. R. Clements, _et al._ "Optimal design for universal multiport
 * interferometers." Optica __3__, pp. 1460-1465 (2016),
 * https://doi.org/10.1364/OPTICA.3.001460
 */
void getClementsDesign(Matrix &M, Point &x, Real eps  = defaultEpsilon);

/** @ingroup design
 * @brief @overload
 *
 * @return Calculated array `x` of phase-shift parameters.
 */
Point getClementsDesign(Matrix &M, Real eps = defaultEpsilon);

} // namespace linopt
