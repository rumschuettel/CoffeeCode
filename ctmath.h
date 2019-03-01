#pragma once

#include <cstddef>
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
	constexpr IntegerT ipow(const IntegerT base, const int exp, const IntegerT result = 1) {
		return exp < 1 ?
			result :
			ipow(base*base, exp / 2, (exp % 2) ? result * static_cast<IntegerT>(base) : result);
	}


	// compile time size of base k tuple
	template<const size_t base, const size_t tuple_length>
	struct BaseKSubsets {
		static constexpr auto count = ipow(base, tuple_length);
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
}