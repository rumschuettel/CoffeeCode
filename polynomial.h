#pragma once

#include "types.h"
#include "ctmath.h"
#include "utility.h"

#include <boost/container_hash/hash.hpp>

#include <map>

#include <assert.h>
// TODO: implement sparse polynomial with larger coefficients

namespace CoffeeCode {
	// monomial types
	// maximum exponent: is equal to the size of the graph,
	//     as we have at most (X wo Y) + (Y wo X) + (X and Y) = (X or Y) <= K_TOT
	//     marked vertices within a set.
	// maximum coefficient: the number of sets giving the same Uidx
	//     can be crudely upper-bounded by 4^K_SYS
	//     This is precisely given by the corresponding MultiplicityT
	struct UnivariateMonomial {
		static constexpr size_t K_TOT = K_SYS + K_ENV;
		static constexpr size_t ExponentWidth = ilog2(K_TOT + 1);

		using ExponentT = SizeStorageType<K_TOT>;
		static constexpr ExponentT MaxExponent = K_TOT;

		using CoefficientT = MultiplicityType<4>;
		using CoefficientArrayT = typename std::array<CoefficientT, MaxExponent + 1>;

		template<typename T>
		inline static ExponentT MakeExponent(const T u1, const T u2, const T u3)
		{
			return checked_cast<ExponentT>(u1 + u2 + u3);
		}

		inline static void AddCoefficientArrays(CoefficientArrayT& lhs, const CoefficientArrayT& rhs)
		{
			for (size_t i = 0; i < lhs.size(); i++)
				lhs[i] += rhs[i];
		}

		inline static void PrintCoefficientArray(std::ostream& stream, const CoefficientArrayT& coefficients)
		{			
			size_t i = coefficients.size();
			for (const auto& coeff : coefficients) {
				stream << +coeff;
				if (--i) stream << ",";
			}
		}

		inline static auto HashCoefficientArray(const CoefficientArrayT& coefficients)
		{
			return boost::hash_range(std::begin(coefficients), std::end(coefficients));
		}
	};

	template<size_t VariableCount = 3>
	struct MultivariateMonomial {
		static constexpr size_t K_TOT = K_SYS + K_ENV;

		using SingleExponentT = SizeStorageType<K_TOT>;
		static constexpr SingleExponentT MaxExponent = K_TOT;
		using ExponentT = std::array<SingleExponentT, VariableCount>;

		using CoefficientT = MultiplicityType<4>;

		// this needs to be an ordered map since we hash the coefficient array
		using CoefficientArrayT = typename std::map<ExponentT, CoefficientT>;

		template<typename T>
		inline static ExponentT MakeExponent(const T& u1, const T& u2, const T& u3)
		{
			return ExponentT{
				checked_cast<SingleExponentT>(u1),
				checked_cast<SingleExponentT>(u2),
				checked_cast<SingleExponentT>(u3)
			};
		}

		inline static void AddCoefficientArrays(CoefficientArrayT& lhs, const CoefficientArrayT& rhs)
		{
			for (const auto& [exponents, coeff] : rhs)
				lhs[exponents] += coeff;
		}

		// print as array elements like [15, [1, 0, 3]]
		// meaning 15 u1^1 u2^0 u3^3
		inline static void PrintCoefficientArray(std::ostream& stream, const CoefficientArrayT& coefficients)
		{			
			size_t i = coefficients.size();
			for (const auto& [exponents, coeff] : coefficients) {
				stream << "[" << +coeff << ",[";
				size_t j = exponents.size();
				for (const auto ex : exponents) {
					stream << +ex;
					if (--j) stream << ",";
				}
				stream << "]]";
				if (--i) stream << ",";
			}
		}

		inline static auto HashCoefficientArray(const CoefficientArrayT& coefficients) noexcept
		{
			size_t seed = 0;
			for (const auto& [exponents, coeff] : coefficients) {
				boost::hash_range(seed, std::begin(exponents), std::end(exponents));
				boost::hash_combine(seed, coeff);
			}
			return seed;
		}
	};

	// polynomial
	template<typename _MonomialT>
	struct Polynomial {
		using MonomialT = _MonomialT;

		using CoefficientT = typename MonomialT::CoefficientT;
		using CoefficientArrayT = typename MonomialT::CoefficientArrayT;
		using ExponentT = typename MonomialT::ExponentT;

		CoefficientArrayT coefficients;

		Polynomial() = default; // zero-initializes coefficients

		void Add(const ExponentT exponent, const CoefficientT coefficient = 1) {
			coefficients[exponent] += coefficient;
		}

		// addition of monomials
		inline Polynomial& operator+=(const ExponentT& rhs) {
			coefficients[rhs] += 1;
			return *this;
		}
		// addition of polynomials
		inline Polynomial& operator+=(const Polynomial& rhs) {
			MonomialT::AddCoefficientArrays(coefficients, rhs.coefficients);
			return *this;
		}

		// to string
		friend std::ostream& operator<< (std::ostream& stream, const Polynomial& poly) {
			MonomialT::PrintCoefficientArray(stream, poly.coefficients);
			return stream;
		}

		// FOR USE IN REDUCE_LAMBDA compressed output
		struct Hash
		{
			inline std::size_t operator()(Polynomial const &poly) const noexcept
			{
				return MonomialT::HashCoefficientArray(poly.coefficients);
			}
		};

		bool operator==(const Polynomial& rhs) const
		{
			return coefficients == rhs.coefficients;
		}
	};

	namespace Std {
	// choose right specialization based on compiler flags
	#ifdef OPTIMIZE_FOR_DEPOLARIZING
		using Polynomial = typename CoffeeCode::Polynomial<UnivariateMonomial>;
	#else
		using Polynomial = typename CoffeeCode::Polynomial<MultivariateMonomial<3>>;
	#endif
	}
}
