#pragma once

#include "ctmath.h"



#include <assert.h>
	// TODO: implement sparse polynomial with larger coefficients

namespace CoffeeCode {
	static constexpr size_t K_TOT = K_SYS + K_ENV;

	// monomial type
	// maximum exponent: is equal to the size of the graph,
	//     as we have at most (X wo Y) + (Y wo X) + (X and Y) = (X or Y) <= K_TOT
	//     marked vertices within a set.
	// maximum coefficient: the number of sets giving the same Uidx
	//     can be upper-bounded by 2^K_TOT, as we have a linear affine map mod 2
	//     K_TOT bits is thus a safe bet.
	struct Monomial {
		static constexpr size_t CoefficientWidth = K_TOT;
		static constexpr size_t ExponentWidth = log2(3*K_TOT);

		using CoefficientT = StdStoreT<CoefficientWidth>;
		using ExponentT = StdStoreT<ExponentWidth>;

		static constexpr CoefficientT MaxCoefficient = 2 << K_TOT;
		static constexpr ExponentT MaxExponent = K_TOT;
	};

	// polynomial
	struct Polynomial {
		using CoefficientT = Monomial::CoefficientT;
		using ExponentT = Monomial::ExponentT;
		static constexpr auto MaxExponent = Monomial::MaxExponent;
		static constexpr auto MaxCoefficient = Monomial::MaxCoefficient;

		CoefficientT coefficients[MaxExponent];

		Polynomial() = default; // zero-initializes coefficients

		void Add(const ExponentT exponent, const CoefficientT coefficient = 1) {
			assert(exponent < MaxExponent);
			coefficients[exponent] += coefficient;
		}

		// addition of polynomials
		inline Polynomial& operator+=(const Polynomial& rhs) {
			for (size_t i = 0; i < MaxExponent; i++)
				coefficients[i] += rhs.coefficients[i];
			return *this;
		}

		// to string
		friend std::ostream& operator<< (std::ostream& stream, const Polynomial& poly) {
			// exploit integer promotion for too-small types
			for (size_t i = 0; i < MaxExponent - 1; i++)
				stream << +poly.coefficients[i] << ", ";
			stream << +poly.coefficients[MaxExponent - 1];
			return stream;
		}
	};
}