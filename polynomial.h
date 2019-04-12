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
		static_assert(VariableCount == 3); // for now only implemented for q1, q2, q3
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

	// numerical monomial with values
	// Vals[1]/Base, Vals[2]/Base, ...
	template<typename FloatT, auto Base, auto... Vals>
	struct NumericalMonomialSinglePoint {
		static constexpr size_t K_TOT = K_SYS + K_ENV;
		static constexpr size_t VariableCount = sizeof...(Vals);
		 // only implemented for 1 or 3 variables
		static_assert(VariableCount == 1 || VariableCount == 3);

		// precompute Vals[1]/Base, ...
		static constexpr std::array<FloatT, VariableCount> FloatVals = {
			(cdiv<FloatT, Vals, Base>) ...
		};

		using SingleExponentT = SizeStorageType<K_TOT>;
		static constexpr SingleExponentT MaxExponent = K_TOT;
		using ExponentT = std::array<SingleExponentT, VariableCount>;

		using CoefficientT = FloatT;

		// encapsulate a value, such that the following proxying takes place:
		// Type[exponent]{value} += coefficient
		// AdditionProxy{exponent, value} += coefficient
		// value += coefficient*FloatVals[exponent]
		using CoefficientArrayT = struct IndexProxy {
			CoefficientT value;
			struct AdditionProxy {
				const ExponentT& exponent;
				CoefficientT& value;
				AdditionProxy(const ExponentT& exponent, CoefficientT& value) 
					: exponent(exponent), value(value)
				{}

				inline void operator+=(const CoefficientT& coeff)
				{
					CoefficientT to_add = coeff;
					for (size_t i = 0; i < VariableCount; i ++)
						to_add *= pow(FloatVals[i], exponent[i]);
					value += to_add;
				}
			};
			inline AdditionProxy operator[](const ExponentT& exponent)
			{
				return AdditionProxy(exponent, value);
			}
			inline bool operator==(const IndexProxy& rhs) const
			{
				return value == rhs.value;
			}
		};

		template<typename T>
		inline static ExponentT MakeExponent(const T u1, const T u2, const T u3)
		{
			if constexpr (VariableCount == 1)
				return ExponentT{ checked_cast<SingleExponentT>(u1 + u2 + u3) };
			else
				return ExponentT{
					checked_cast<SingleExponentT>(u1),
					checked_cast<SingleExponentT>(u2),
					checked_cast<SingleExponentT>(u3)
				};
		}

		inline static void AddCoefficientArrays(CoefficientArrayT& lhs, const CoefficientArrayT& rhs)
		{
			lhs.value += rhs.value;
		}

		inline static void PrintCoefficientArray(std::ostream& stream, const CoefficientArrayT& value)
		{			
			stream << value.value;
		}

		// hashing float values is generally a bad idea
		// since this is just used to compress the output,
		// we don't really care.
		inline static size_t HashCoefficientArray(const CoefficientArrayT& value)
		{
			// no data race if only using the const interface
			static std::hash<FloatT> hasher;
			return hasher(value.value);
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

		inline bool operator==(const Polynomial& rhs) const
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
