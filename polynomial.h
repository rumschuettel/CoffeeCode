#pragma once

#include "ctmath.h"

#include <boost/container_hash/hash.hpp>

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
	//     K_TOT bits is thus a safe bet; since we actually might need to store 2^K_TOT, we add 1
	struct Monomial {
		static constexpr size_t CoefficientWidth = K_TOT + 1;
		static constexpr size_t ExponentWidth = ilog2(K_TOT + 1);

		using CoefficientT = StdStoreT<CoefficientWidth>;
		using ExponentT = StdStoreT<ExponentWidth>;

		static constexpr CoefficientT MaxCoefficient = Bitmask<CoefficientT, K_TOT>::mask1000;
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

		// hash
		struct Hash
		{
			inline std::size_t operator()(Polynomial const &poly) const noexcept
			{
				return boost::hash_range(
					std::begin(poly.coefficients),
					std::end(poly.coefficients)
				);
			}
		};

		// comparison for this type
		bool operator==(const Polynomial& rhs) const
		{
			return std::equal(
				std::begin(coefficients),
				std::end(coefficients),
				std::begin(rhs.coefficients)
			);
		}
	};
}
