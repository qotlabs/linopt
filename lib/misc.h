// SPDX-License-Identifier: LGPL-3.0-or-later
// SPDX-FileCopyrightText: Copyright (c) 2018-2025, Quantum Optical Technologies Laboratories
// SPDX-FileContributor: Struchalin Gleb <struchalin.gleb@physics.msu.ru>
// SPDX-FileContributor: Dyakonov Ivan <iv.dyakonov@physics.msu.ru>

#pragma once

#include "types.h"
#include <cmath>
#include <ostream>

namespace linopt
{

#ifdef __GNUC__
	static inline int ctz(unsigned x)
	{
		return __builtin_ctz(x);
	}
	static inline int ctz(unsigned long x)
	{
		return __builtin_ctzl(x);
	}
	static inline int ctz(unsigned long long x)
	{
		return __builtin_ctzll(x);
	}
#elif _MSC_VER
	#include <intrin.h>
	#pragma intrinsic(_BitScanForward)
	static inline int ctz(unsigned long x)
	{
		unsigned long i;
		_BitScanForward(&i, x);
		return i;
	}
	static inline int ctz(__int64 x)
	{
		unsigned long i;
		_BitScanForward64(&i, x);
		return i;
	}
#else
	#pragma message("Builtin ctz() is not available. Using C version of ctz().")
	static inline int ctz(uint16_t x)
	{
		int n = 1;
		if ((x & 0x00FF) == 0) {n += 8;	x >>= 8;}
		if ((x & 0x000F) == 0) {n += 4;	x >>= 4;}
		if ((x & 0x0003) == 0) {n += 2;	x >>= 2;}
		return n - (x & 1);
	}
	static inline int ctz(uint32_t x)
	{
		int n = 1;
		if ((x & 0x0000FFFF) == 0) {n += 16; x >>= 16;}
		if ((x & 0x000000FF) == 0) {n += 8;  x >>= 8;}
		if ((x & 0x0000000F) == 0) {n += 4;  x >>= 4;}
		if ((x & 0x00000003) == 0) {n += 2;  x >>= 2;}
		return n - (x & 1);
	}
	static inline int ctz(uint64_t x)
	{
		int n = 1;
		if ((x & 0x00000000FFFFFFFF) == 0) {n += 32; x >>= 32;}
		if ((x & 0x000000000000FFFF) == 0) {n += 16; x >>= 16;}
		if ((x & 0x00000000000000FF) == 0) {n += 8;  x >>= 8;}
		if ((x & 0x000000000000000F) == 0) {n += 4;  x >>= 4;}
		if ((x & 0x0000000000000003) == 0) {n += 2;  x >>= 2;}
		return n - (x & 1);
	}
#endif

template<typename T>
inline T mod(T x, T a)
{
	x = std::fmod(x, a);
	return (x < 0) ? x + a : x;
}

inline int isqrt(int x)
{
	return static_cast<int>(std::sqrt(x) + 0.5);
}

template<typename T>
std::ostream& printArray(std::ostream &stream, const T &a,
						 const char *b1 = "{",
						 const char *delim = ", ",
						 const char *b2 = "}")
{
	if(a.empty())
	{
		stream << b1 << b2;
	}
	else
	{
		stream << b1;
		for(auto iter = a.begin(); iter != --a.end(); iter++)
			stream << *iter << delim;
		stream << *(--a.end()) << b2;
	}
	return stream;
}

std::ostream& printComplex(std::ostream &stream, const Complex &x);

} // namespace linopt
