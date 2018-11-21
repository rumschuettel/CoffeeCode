#pragma once

// monomial
struct Monomial {
	using CoefficientT = uint32_t;
	using ExponentT = uint8_t;
	static constexpr auto MaxCoefficient = std::numeric_limits<CoefficientT>::max();
};

// polynomial
template<const size_t max_exponent>
struct Polynomial {
	Monomial::CoefficientT coefficients[max_exponent];

	Polynomial() = default; // zero-initializes coefficients

	void Add(const Monomial::ExponentT exponent) {
		assert(exponent < max_exponent);
		coefficients[exponent] ++;
	}
	void Add(const Polynomial& other) {
		for (size_t i = 0; i < max_exponent; i++) {
			assert(coefficients[i] + other.coefficients[i] < Monomial::MaxCoefficient);
			coefficients[i] += other.coefficients[i];
		}
	}
};
