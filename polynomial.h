#pragma once

#include <assert.h>
// TODO: implement sparse polynomial with larger coefficients

// monomial
template<typename _CoefficientT = uint32_t, typename _ExponentT = uint8_t>
struct Monomial {
	using CoefficientT = _CoefficientT;
	using ExponentT = _ExponentT;
	static constexpr auto MaxCoefficient = std::numeric_limits<CoefficientT>::max();
	static constexpr auto MaxExponent = std::numeric_limits<ExponentT>::max();
};

// polynomial
template<const size_t max_exponent>
struct Polynomial {
	Monomial<>::CoefficientT coefficients[max_exponent];

	Polynomial() = default; // zero-initializes coefficients

	void Add(const Monomial<>::ExponentT exponent, const Monomial<>::CoefficientT coefficient = 1) {
		assert(exponent < max_exponent);
		coefficients[exponent] += coefficient;
	}
	void Add(const Polynomial& other) {
		for (size_t i = 0; i < max_exponent; i++) {
			assert(coefficients[i] + other.coefficients[i] < Monomial<>::MaxCoefficient);
			coefficients[i] += other.coefficients[i];
		}
	}

	// to string
	friend std::ostream& operator<< (std::ostream& stream, const Polynomial<max_exponent>& poly) {
		for (size_t i = 0; i < max_exponent - 1; i++)
			stream << poly.coefficients[i] << ", ";
		stream << poly.coefficients[max_exponent - 1];
		return stream;
	}
};
