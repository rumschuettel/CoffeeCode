#include "CoffeeCode.h"

#include <iostream>
#include <cstdint>
#include <vector>
#include <assert.h>

using std::cout;
using std::endl;

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
template<const size_t number_of_1s>
struct Bitmask1s {
	static constexpr auto mask = (1ull << number_of_1s) - 1;
};

// monomial
struct Monomial {
	using CoefficientT = uint16_t;
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


int main()
{
	constexpr auto k_sys = 18, k_env = 1, k_tot = k_sys + k_env;
	CoffeeCode::AdjacencyMatrix<k_sys, k_env> M{ "0111111111110000000100000000000000000010000000000000000001000000000000\
0000001000000000000000000100000000000000000010000000000000000001000000\
0000000000001000000000000000000100000000000000000010000000000000000001\
0000000000011111110000000000010000000000000000001000000000000000000100\
0000000000000000100000000000000000010000000000000000001000000000000000\
00010000000" };

	// Calculate lambda and lambda_pre
	constexpr auto max_exponent = k_tot;
	std::vector<Polynomial<max_exponent>> lambda(ipow(2, k_tot)), lambda_pre(ipow(2, k_sys));

	// iterate over all U1, U2 and U3 subsets
	// we do not cache as memory access is potentially slower than counting
	using VectorT = decltype(M)::RowVectorT;
	using SubsetT = VectorT::StoreT;

	// could reach numerical limits, hence break condition at end

	for (SubsetT subsetX = 0; ; subsetX++) {
		for (SubsetT subsetY = 0; ; subsetY++) {
			// we explicitly allow downcasting, so no {} as initializer; in fact the compiler complains if we do.
			// note that the types match, though, as subsetX is of type VectorT::StoreT, and the bit operations
			// keep it within the allowed range.
			const auto XwoY = VectorT(subsetX & ~subsetY);
			const auto YwoX = VectorT(subsetY & ~subsetX);
			const auto XnY = VectorT(subsetX & subsetY);

			const Monomial::ExponentT u1 = static_cast<Monomial::ExponentT>(XwoY.popcount());
			const Monomial::ExponentT u2 = static_cast<Monomial::ExponentT>(YwoX.popcount());
			const Monomial::ExponentT u3 = static_cast<Monomial::ExponentT>(XnY.popcount());

			const SubsetT Uidx = (M * (XwoY + YwoX) + (YwoX + XnY)).vec;
			const SubsetT UAidx = Uidx & Bitmask1s<k_sys>::mask;

			// add monomials
			assert(u1 + u2 + u3 < max_exponent);
			const Monomial::ExponentT p_exponent = u1 + u2 + u3;
  			lambda[Uidx].Add(p_exponent);
			lambda_pre[UAidx].Add(p_exponent);

			if (subsetY == BaseKSubsets<2, k_sys>::count-1) break;
		}

		if (subsetX == BaseKSubsets<2, k_sys>::count-1) break;
	}

	// Calculate lambda_a
	constexpr auto max_exponent_summed = k_tot;
	std::vector<Polynomial<max_exponent_summed>> lambda_a(ipow(2, k_sys));

	// build Ulookup
	const auto MAB = M.AB();
	using SubsetAT = decltype(MAB)::ColumnVectorT::StoreT;
	using SubsetBT = decltype(MAB)::RowVectorT::StoreT;
	for (SubsetAT subsetA = 0; ; subsetA++) {
		for (SubsetBT subsetB = 0; ; subsetB++) {
			const auto Ulookup = (MAB * subsetB + subsetA).vec;
			lambda_a[subsetA].Add(lambda_pre[Ulookup]);

			if (subsetB == BaseKSubsets<2, k_env>::count-1) break;
		}

		if (subsetA == BaseKSubsets<2, k_sys>::count-1) break;
	}

	// dump output to Mathematica-readable file
	cout << "{{";
	for (const auto& polynomial : lambda) {
		cout << "{";
		for (const auto& coefficient : polynomial.coefficients) {
			cout << static_cast<size_t>(coefficient) << ", ";
		}
		cout << "Nothing},";
	}
	cout << "Nothing},\n{";
	for (const auto& polynomial : lambda_a) {
		cout << "{";
		for (const auto& coefficient : polynomial.coefficients) {
			cout << static_cast<size_t>(coefficient) << ", ";
		}
		cout << "Nothing},";
	}
	cout << "Nothing}}" << endl;

	return 0;
}
