#include "CoffeeCode.h"

#include <unordered_map>
#include <unordered_set>
#include <chrono>
#include <bitset>

namespace CoffeeCode {
	// extract compile time parameters
	template<typename T>
	struct SymmetricInstance {
		// parameters given
		constexpr static size_t k_sys = K_SYS;
		constexpr static size_t k_env = K_ENV;
		constexpr static size_t k_tot = k_sys + k_env;
		using sgs = typename T::sgs;

		// validate parameters
		static_assert(T::sgs::Length == k_sys);
		static_assert(T::adjacency_matrix.size() == k_tot);
		static_assert(T::adjacency_matrix[0].size() == k_tot);

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

	#include "cc-instance.h"

	using instance = CoffeeCode::SymmetricInstance<graphstate_instance>;
}


// nauty library
#include "nautylink.h"

int SymmetricSolver() {
	using VectorT = typename instance::RowVectorT;
	using SubsetT = typename VectorT::StoreT;

	using ExponentT = typename CoffeeCode::Monomial::ExponentT;
	using CoefficientT = typename CoffeeCode::Monomial::CoefficientT;

	auto nauty = CoffeeCode::NautyLink::NautyLink(instance::M);
	using CanonicalImageT = decltype(nauty)::CanonicalImageT;

	const auto fullGroupOrder = nauty.GroupOrder();

	auto start = std::chrono::steady_clock::now();

	// CHANNEL ACTION

	// lambdas as hash maps
	// TODO: replace with https://github.com/greg7mdp/sparsepp
	using LambdaT = std::unordered_map<
		CanonicalImageT,
		std::tuple<
			CoffeeCode::Polynomial,
			CoffeeCode::NautyLink::OrbitSizeT
		>,
		CanonicalImageT::Hash
	>;
	LambdaT lambda, lambda_pre;

	size_t counter = 0;
	for (const auto& tuple : instance::sgs::TupleCosets<4>()) {
		counter++;

		// calculate multiplicity of base 4 tuple
		nauty.SetColoring(tuple);
		const auto stabGroupOrder4 = nauty.GroupOrder();
		const auto orbitSize4 = fullGroupOrder / stabGroupOrder4;

		// get low and high bit from tuples
		// note that SubsetT is such that the i'th index equals a bit shit i to the left
		SubsetT subsetX{ 0 }, subsetY{ 0 };
		for (size_t i = 0; i < instance::k_sys; i++) {
			const auto low_bit = tuple[i] & 0b01;
			const auto high_bit = (tuple[i] & 0b10) >> 1;

			CoffeeCode::OrBit(subsetX, !!low_bit, i);
			CoffeeCode::OrBit(subsetY, !!high_bit, i);
		}

		// apply channel
		const auto term = ChannelAction(subsetX, subsetY, instance::M);
		
		// monomial exponents
		const ExponentT p_exponent = term.uSum();

		//// A: Add to lambda
		{
			// multiplicity of resulting base 2 tuple
			const auto[UIdxCanonical, stabGroupOrder2] = nauty.CanonicalColoring(term.Uidx);
			const auto orbitSize2 = fullGroupOrder / stabGroupOrder2;
			const CoefficientT coeff = static_cast<CoefficientT>(orbitSize4 / orbitSize2);

			// accumulate polynomial with potentially pre-existing terms
			// note that mult is equal for U123 that map to the same UIdxCanonical
			auto& [poly, mult] = lambda[UIdxCanonical];
			poly.Add(p_exponent, coeff);
			mult = orbitSize2;
		}

		//// B: Add to lambda_pre
		{
			// we project out the environment for term.Uidx
			const auto[UIdxCanonical_pre, stabGroupOrder2_pre] = nauty.CanonicalColoring(
				term.Uidx & CoffeeCode::Bitmask<decltype(term.Uidx), instance::k_sys>::mask0111
			);
			const auto orbitSize2_pre = fullGroupOrder / stabGroupOrder2_pre;
			const CoefficientT coeff = static_cast<CoefficientT>(orbitSize4 / orbitSize2_pre);

			// accumulate polynomial with same terms as above
			// multiplicity is 1
			auto& [poly, mult] = lambda_pre[UIdxCanonical_pre];
			poly.Add(p_exponent, coeff);
			mult = orbitSize2_pre;
		}
	}


	// PARTIAL TRACE
	using SubsetBT = decltype(instance::MAB)::RowVectorT::StoreT;
	using SubsetAT = decltype(instance::MAB)::ColumnVectorT::StoreT;
	LambdaT lambda_a;

	for (const auto& tuple : instance::sgs::TupleCosets<2>()) {
		// TupleT to SubsetT and HashT because that one can be different
		SubsetAT subsetA{ 0 };
		for (size_t i = 0; i < instance::k_sys; i++) 
			CoffeeCode::OrBit(subsetA, !!tuple[i], i);
		const auto [Akey, _] = nauty.CanonicalColoring(subsetA);

		// this inner loop is as in CCFull with the break condition at the end to avoid overflows
		for (SubsetBT subsetB = 0; ; subsetB++) {
			// new key, subset in A but cast to full subset on entire graph
			const auto BtoA = static_cast<SubsetT>((instance::MAB * subsetB + subsetA).vec);

			// find canonical image for key
			const auto [key, keyStabOrder] = nauty.CanonicalColoring(BtoA);
			const auto keyMult = fullGroupOrder / keyStabOrder;

			//// C: Add to lambda_a
			const auto& [poly_pre, mult_pre] = lambda_pre[key];
			auto& [poly_a, mult_a] = lambda_a[Akey];

			assert(keyMult == mult_pre);  // must hold by definition
			poly_a += poly_pre;
			mult_a = keyMult;
			
			if (subsetB == CoffeeCode::BaseKSubsets<2, instance::k_env>::count - 1) break;
		}
	}


	// EXPORT AS JSON
	// some statistics
	std::cout << "{\n\"tuples\": " << counter << ",\n\"time\": ";
	auto end = std::chrono::steady_clock::now();
	std::cout << std::chrono::duration <double, std::milli> (end-start).count() << ",\n";

	// lambda
	auto print_lambda = [&](auto lambda) -> void {
		size_t i = lambda.size();
		for (const auto& [key, value] : lambda) {
			const auto& [poly, mult] = value;
			std::cout << "[" << mult << ", [" << poly << "]]";
			if (--i) std::cout << ",";
			std::cout << "\n";
		}
	};
	std::cout << "\"lambda\": [\n";
	print_lambda(lambda);
	std::wcout << "],\n";
	// lambda_a
	std::cout << "\"lambda_a\": [\n";
	print_lambda(lambda_a);
	std::cout << "]\n}\n";

	return RET_OK;
}