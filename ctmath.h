#pragma once

#include "ctlookup.h"
#include "utility.h"

#include <cstddef>
#include <stdint.h>

#include <limits>
#include <cmath>

// we need to define the wider types here in case we specialize any math function
#include <boost/multiprecision/cpp_int.hpp>

using boost::multiprecision::uint128_t;
using boost::multiprecision::uint256_t;
using boost::multiprecision::uint512_t;

// https://github.com/calccrypto/uint256_t.git
// provides fast 128 and 256 bit integers
//#include "src/vectorclass/vectorclass.h"


namespace LibPopcount {
	#include "libpopcnt.h"
}

namespace CoffeeCode {
	// bit fiddling
	inline size_t Popcount(const uint64_t& s) noexcept
	{
		return LibPopcount::popcount64(s);
	}
	template<typename T>
	inline size_t Popcount(const boost::multiprecision::number<T>& s)
	{
		const auto* limbs = s.backend().limbs();
		const size_t limb_count = s.backend().size();
		return LibPopcount::popcnt(limbs, limb_count * sizeof(limbs[0]));
	}

	template<typename StoreT>
	constexpr inline void OrBit(StoreT& s, const bool bit, const size_t i)
	{
		s |= static_cast<StoreT>(bit) << i;
	}

	// compile-time ceil(log2)
    constexpr size_t ilog2(const size_t n) {
        return n <= 1 ? 0 : 1 + ilog2((n + 1) / 2);
    }

	// approximate ceil(log2(n!))
	constexpr size_t ilog2factorial(const size_t n) {
		return n <= 1 ? 1 : ((n + 1)*ilog2(n) - n);
    }

	// compile time integer exponent
	template<typename IntegerT>
	constexpr IntegerT ipow(const IntegerT base, const IntegerT exp, const IntegerT result = 1) {
		return exp < 1 ?
			result :
			ipow(base*base, exp / 2, (exp % 2) ? result*base : result);
	}


	// compile time size of base k tuple
	template<const size_t base, const size_t tuple_length>
	struct BaseKSubsets {
		static constexpr size_t count = ipow(base, tuple_length);
	};

	// compile time bitmask with k 1s
	template<
		typename MaskT,
		size_t number_of_1s,
		typename = std::enable_if_t<std::numeric_limits<MaskT>::digits >= number_of_1s + 1>
	>
	struct Bitmask {
		static constexpr auto mask1000 = LUTs::Bitmasks<MaskT>::lut1000[number_of_1s];
		static constexpr auto mask0111 = LUTs::Bitmasks<MaskT>::lut0111[number_of_1s];
	};

	// n choose k
	template<typename SizeT>
	constexpr inline SizeT nCHk(size_t n, size_t k) {
		if (k == 0) return 1;
		return n * nCHk<SizeT>(n - 1, k - 1) / k;
	}


	// factorial function using lookup tables
	template<typename SizeT>
	constexpr inline SizeT Factorial(size_t n)
	{
		assert(n < sizeof(LUTs::Factorial<SizeT>::lut));
		return LUTs::Factorial<SizeT>::lut[n];
	}

	// binomial function using lookup tables
	// precalculate on startup
	template<typename SizeT>
	inline SizeT Binomial(const size_t n, const size_t k)
	{
		constexpr size_t LUTSize = std::numeric_limits<SizeT>::digits;
		static SizeT lut[LUTSize][LUTSize];
		static bool filled = false;

		if (!filled) {
			for (size_t n = 0; n < LUTSize; n++)
				for (size_t k = 0; k <= n; k++)
					lut[n][k] = nCHk<SizeT>(n, k);

			filled = true;
		}

		assert(n < LUTSize);
		assert(k <= n);

		return lut[n][k];
	}
	
	// or use built-in lookup table for small types
#define BINOMIAL_LUT(TYPE)\
	template<>\
	inline TYPE Binomial<TYPE>(const size_t n, const size_t k) noexcept\
	{\
		constexpr auto& lut = LUTs::BinomialCoefficient<TYPE>::lut;\
\
		assert(n < sizeof(lut));\
		assert(k < sizeof(lut[n]));\
		assert(k <= n);\
\
		return lut[n][k];\
	}

	BINOMIAL_LUT(uint8_t)
	BINOMIAL_LUT(uint16_t)
	BINOMIAL_LUT(uint32_t)
	BINOMIAL_LUT(uint64_t)
#undef BINOMIAL_LUT

	template<>
	inline double Binomial<double>(const size_t n, const size_t k)
	{
		constexpr auto& lut = LUTs::BinomialCoefficient<double>::lut;
		
		assert(k <= n);
		if (n < sizeof(lut)) {
			assert(k < sizeof(lut[n]));
			return lut[n][k];
		}

		// otherwise approximate
		return boost::math::binomial_coefficient<double>(checked_cast<unsigned int>(n), checked_cast<unsigned int>(k));
	}
}