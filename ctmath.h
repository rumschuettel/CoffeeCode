#pragma once

#include "ctlookup.h"

#include <cstddef>
#include <stdint.h>

#include <limits>
#include <cmath>

namespace CoffeeCode {
	// compile-time ceil(log2)
    constexpr size_t ilog2(const size_t n) {
        return n <= 1 ? 0 : 1 + ilog2((n + 1) / 2);
    }

	// approximate ceil(log2(n!))
	constexpr size_t ilog2factorial(const size_t n) {
		return (n + 1)*ilog2(n) - n;
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
        static constexpr auto mask1000 = static_cast<MaskT>( 0b01 ) << number_of_1s;
		static constexpr auto mask0111 = mask1000 - 1;
	};


	// factorial function using lookup tables
	template<typename SizeT>
	constexpr inline SizeT Factorial(size_t n)
	{
		assert(n < sizeof(LUTs::Factorial<SizeT>::lut));
		return LUTs::Factorial<SizeT>::lut[n];
	}

	// binomial function using lookup tables
	template<typename SizeT>
	constexpr inline SizeT Binomial(const SizeT n, const SizeT k)
	{
		assert(n < sizeof(LUTs::BinomialCoefficient<SizeT>::lut));
		assert(k < sizeof(LUTs::BinomialCoefficient<SizeT>::lut[n]));
		return LUTs::BinomialCoefficient<SizeT>::lut[n][k];
	}
}