#pragma once

#include "ctmath.h"



#include <assert.h>
	// TODO: implement sparse polynomial with larger coefficients

namespace CoffeeCode {
	// monomial type
	// maximum exponent is 3 * the size of the graph, so uint16_t should be plenty of room
	// maximum coefficient: the number of sets giving the same Uidx
	//     can be upper-bounded by the number of sets overall, which is 4^{K_SYS + K_ENV}.
	//     2 * K_TOT bits is thus a safe bet.
	struct Monomial {
		static constexpr size_t CoefficientWidth = 2 * (K_SYS + K_ENV);
		static constexpr size_t ExponentWidth = log2(3 * (K_SYS + K_ENV));

		using CoefficientT = StdStoreT<CoefficientWidth>;
		using ExponentT = StdStoreT<ExponentWidth>;
		static constexpr auto MaxCoefficient = std::numeric_limits<CoefficientT>::max();
		static constexpr auto MaxExponent = std::numeric_limits<ExponentT>::max();
	};

	// polynomial
	template<const size_t max_exponent>
	struct Polynomial {
		Monomial::CoefficientT coefficients[max_exponent];

		Polynomial() = default; // zero-initializes coefficients

		void Add(const Monomial::ExponentT exponent, const Monomial::CoefficientT coefficient = 1) {
			assert(exponent < max_exponent);
			coefficients[exponent] += coefficient;
		}
		void Add(const Polynomial& other) {
			for (size_t i = 0; i < max_exponent; i++) {
				assert(coefficients[i] + other.coefficients[i] < Monomial::MaxCoefficient);
				coefficients[i] += other.coefficients[i];
			}
		}

		// to string
		friend std::ostream& operator<< (std::ostream& stream, const Polynomial<max_exponent>& poly) {
			// exploit integer promotion for too-small types
			for (size_t i = 0; i < max_exponent - 1; i++)
				stream << +poly.coefficients[i] << ", ";
			stream << +poly.coefficients[max_exponent - 1];
			return stream;
		}
	};
}