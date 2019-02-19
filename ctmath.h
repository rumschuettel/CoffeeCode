#pragma once

#include <cstddef>
#include <limits>
#include <cmath>

namespace CoffeeCode {
	// compile-time ceil(log2)
    constexpr size_t log2(const size_t n) {
        return n <= 1 ? 0 : 1 + log2((n + 1) / 2);
    }
    // compile-time log-factorial, i.e. lf(n) = log_2(n!) = log_2(n) + log_2(n-1) + ... + log_2(1)
    constexpr double log2factorial(const double n) {
        return n <= 1 ? 0 : ::log2(n) + log2factorial(n-1);
    }
    template<size_t n>
    constexpr size_t ilog2factorial() {
        constexpr double dn = static_cast<double>(n);
        constexpr double l2f = log2factorial(dn);
        return static_cast<size_t>(::ceil(l2f));
    }

	// compile time integer exponent
	constexpr size_t ipow(const size_t base, const int exp, const size_t result = 1) {
		return exp < 1 ? result : ipow(base*base, exp / 2, (exp % 2) ? result * base : result);
	}


	// compile time size of base k tuple
	template<const size_t base, const size_t tuple_length>
	struct BaseKSubsets {
		static constexpr auto count = ipow(base, tuple_length);
	};
	// compile time bitmask with k 1s
	template<typename MaskT, const size_t number_of_1s>
	struct Bitmask {
        static constexpr auto mask1000 = static_cast<MaskT>( 0b01 ) << number_of_1s;
		static constexpr auto mask0111 = mask1000 - 1;
	};
}