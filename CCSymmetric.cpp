#include "CoffeeCode.h"

#include <unordered_map>
#include <chrono>
#include <bitset>
#include <map>
#include <set>

namespace CoffeeCode {
	// extract compile time parameters
	template<typename T>
	struct SymmetricInstance {
		// parameters given
		constexpr static size_t k_sys = T::k_sys;
		constexpr static size_t k_env = T::k_env;
		constexpr static size_t k_tot = k_sys + k_env;
		using sgs = typename T::sgs;
		static_assert(T::sgs::Length == k_sys);

		// adjacency matrix
		using MatrixT = CoffeeCode::AdjacencyMatrix<k_sys, k_env>;
		constexpr static auto M = MatrixT{ T::adjacency_matrix };
		constexpr static auto MAB = M.AB();

		// types
		using RowVectorT = typename MatrixT::RowVectorT;
	};

}

namespace {
	enum ReturnValue {
		RET_OK = 0,
		RET_WRONG_INPUT = 1
	};


	using CoffeeCode::SGSTransversal;
	using CoffeeCode::SGSGenerator;
	using CoffeeCode::Group;
	using CoffeeCode::Permutation;
	using CoffeeCode::AdjacencyMatrix;

	template<size_t S>
	using AdjacencyMatrixT = std::array<std::array<CoffeeCode::StdBitT, S>, S>;

	// this is a compile-time problem instance
	struct graphstate_instance {
		using sgs = SGSTransversal<
			SGSGenerator<1, Group<
			Permutation<2, 1, 0>
			>>
			>;
		constexpr static size_t k_sys = 3, k_env = 1;
		constexpr static AdjacencyMatrixT<4> adjacency_matrix{ {{0, 1, 0, 0}, {1, 0, 1, 1}, {0, 1, 0, 0}, {0, 1, 0, 0}} };
	};


	using instance = CoffeeCode::SymmetricInstance<graphstate_instance>;
}


// nauty library
#include "nautylink.h"

int SymmetricSolver() {
	using VectorT = typename instance::RowVectorT;
	using SubsetT = typename VectorT::StoreT;

	using ExponentT = typename Monomial<>::ExponentT;
	using CoefficientT = typename Monomial<>::CoefficientT;

	constexpr auto max_exponent = instance::k_tot;

	auto nauty = CoffeeCode::NautyLink::NautyLink(instance::M);
	const auto fullGroupOrder = nauty.GroupOrder();

	auto start = std::chrono::steady_clock::now();

	// CHANNEL ACTION

	// lambdas as hash maps
	// TODO: replace with https://github.com/greg7mdp/sparsepp
	using LambdaT = std::unordered_map<
		decltype(nauty)::CanonicalImageT,
		std::tuple<Polynomial<max_exponent>, size_t, SubsetT>,
		decltype(nauty)::CanonicalImageT::Hash
	>;
	LambdaT lambda;

	std::map<size_t, std::set<size_t>> debug;

	size_t counter = 0;
	for (const auto& tuple : instance::sgs::TupleCosets<4>()) {
		counter++;

		// calculate multiplicity of base 4 tuple
		nauty.SetColoring(tuple);
		const auto orbitSize4 = fullGroupOrder / nauty.GroupOrder();

		// get low and high bit from tuples
		SubsetT subsetX{ 0 }, subsetY{ 0 };
		for (size_t i = 0; i < instance::k_sys; i++) {
			subsetX |= static_cast<SubsetT>(tuple[i] & 0b01) << i;
			subsetY |= static_cast<SubsetT>(tuple[i] & 0b10) << (i - 1); // -1 because stored one bit up anyhow
		}
		// pad environment vertices on right
		subsetX <<= instance::k_env;
		subsetY <<= instance::k_env;

		// apply channel and get canonical image
		const auto term = ChannelAction(subsetX, subsetY, instance::M);
		const auto[UidxCanonical, stabGroupOrder2] = nauty.CanonicalColoring(term.Uidx);
		// multiplicity of resulting base 2 tuple
 		const auto orbitSize2 = fullGroupOrder / stabGroupOrder2;

		// add monomials
		const ExponentT p_exponent = term.u1 + term.u2 + term.u3;
		assert(p_exponent < max_exponent);
		const CoefficientT coeff = static_cast<CoefficientT>(orbitSize4 / orbitSize2);
		auto& [poly, mult, original_Uidx] = lambda[UidxCanonical];
		poly.Add(p_exponent, coeff);
		mult = orbitSize2; // mult is, by construction, always identical for identical UidxCanonical
		original_Uidx = term.Uidx; // we store one representative, it doesn't matter which one

		print(tuple);
		auto hashF = decltype(nauty)::CanonicalImageT::Hash{};
		std::cout << std::bitset<4>(term.Uidx) << " " << (hashF(UidxCanonical)%1000) << "\n";

		debug[hashF(UidxCanonical)].insert(term.Uidx);
	}

	for (auto idx : debug) {
		std::cout << (idx.first % 1000) << "\n";
		for (auto el : idx.second)
			std::cout << std::bitset<4>(el) << "\n";
		std::cout << "\n";
	}

	// PARTIAL TRACE
	using SubsetBT = decltype(instance::MAB)::RowVectorT::StoreT;
	using SubsetAT = decltype(instance::MAB)::ColumnVectorT::StoreT;
	using LambdaAT = std::unordered_map<
		SubsetT,
		std::tuple<Polynomial<max_exponent>, size_t>
	>;
	LambdaAT lambda_a;

	for (const auto& l : lambda) {
		const auto& [poly, mult, original_Uidx] = l.second;
		// extract system bits
		const auto subsetA = static_cast<SubsetAT>((original_Uidx >> instance::k_env) & Bitmask1s<instance::k_sys>::mask);

		// this inner loop is as in CCFull with the break condition at the end to avoid overflows
		for (SubsetBT subsetB = 0; ; subsetB++) {
			const auto Ulookup = (instance::MAB * subsetB + subsetA).vec;

			// TODO: canonicalize Ulookup
			auto& [poly_a, mult_a] = lambda_a[Ulookup];
			poly_a.Add(poly);
			mult_a = mult;

			if (subsetB == BaseKSubsets<2, instance::k_env>::count - 1) break;
		}
	}


	auto end = std::chrono::steady_clock::now();
	std::cout << "\ncounted " << counter << " tuples.\n";
	std::cout << std::chrono::duration <double, std::milli> (end-start).count() << " ms\n\n";

	return RET_OK;
}