#pragma once

#include "types.h"
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
	//     can be upper-bounded by 2^K_SYS, as we have a linear affine map mod 2
	//     starting from {0, 1, 2, 3}^K_SYS.
	//     This is precisely given by the corresponding MultiplicityT
	struct Monomial {
		static constexpr size_t ExponentWidth = ilog2(K_TOT + 1);

		using CoefficientT = MultiplicityType<4>;
		using ExponentT = StdStoreT<ExponentWidth>;

		static constexpr ExponentT MaxExponent = K_TOT;
	};

	// polynomial
	struct Polynomial {
		using CoefficientT = Monomial::CoefficientT;
		using ExponentT = Monomial::ExponentT;
		static constexpr auto MaxExponent = Monomial::MaxExponent;

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
