#include "CoffeeCode.h"

#include <unordered_map>
#include <chrono>

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
	using InitializerArrayT = std::array<std::array<CoffeeCode::StdBitT, S>, S>;

	// this is a compile-time problem instance
	struct double_legged_graph {
		using sgs = SGSTransversal<
			SGSGenerator<2, Group<
			Permutation<0, 4, 5, 6, 1, 2, 3>
			>>,
			SGSGenerator<3, Group<
			Permutation<0, 1, 3, 2, 4, 6, 5>
			>>,
			SGSGenerator<6, Group<
			Permutation<0, 1, 2, 3, 4, 6, 5>
			>>
			>;

		constexpr static size_t k_sys = 7, k_env = 1;
		constexpr static InitializerArrayT<8> adjacency_matrix{ { {0, 1, 0, 0, 1, 0, 0, 1}, {1, 0, 1, 1, 0, 0, 0, 0}, {0, 1, 0, 0, 0, 0, 0, 0}, {0, 1, 0, 0, 0, 0, 0, 0}, {1, 0, 0, 0, 0, 1, 1, 0}, {0, 0, 0, 0, 1, 0, 0, 0}, {0, 0, 0, 0, 1, 0, 0, 0}, {1, 0, 0, 0, 0, 0, 0, 0} } };
	};


	using instance = CoffeeCode::SymmetricInstance<double_legged_graph>;
}


// nauty library
#include "nautylink.h"

int SymmetricSolver() {
	using VectorT = typename instance::RowVectorT;
	using SubsetT = typename VectorT::StoreT;

	using ExponentT = typename Monomial<>::ExponentT;

	constexpr auto max_exponent = instance::k_tot;

	// lambdas as hash maps
	// TODO: replace with https://github.com/greg7mdp/sparsepp
	using LambdaT = std::unordered_map<SubsetT, Polynomial<max_exponent>>;
	LambdaT lambda;
	LambdaT lambda_a;

	auto nauty = CoffeeCode::NautyLink::NautyLink(instance::M);
	std::cout << nauty.GroupOrder();

	return RET_OK;

	auto start = std::chrono::steady_clock::now();

	size_t counter = 0;
	for (const auto& tuple : instance::sgs::TupleCosets<4>()) {
		if (counter++ % 1000 == 0)
			print(tuple);
		
		continue;
		// get low and high bit from tuples
		SubsetT subsetX{ 0 }, subsetY{ 0 };
		for (size_t i = 0; i < instance::k_sys; i++) {
			subsetX |= static_cast<SubsetT>(tuple[i] & 0b01) << i;
			subsetY |= static_cast<SubsetT>(tuple[i] & 0b10) << (i - 1); // -1 because stored one bit up anyhow
		}

		// apply channel
		const auto term = ChannelAction(subsetX, subsetY, instance::M);
		const SubsetT UAidx = term.Uidx & Bitmask1s<instance::k_sys>::mask;

		// add monomials
		assert(term.u1 + term.u2 + term.u3 < max_exponent);
		const ExponentT p_exponent = term.u1 + term.u2 + term.u3;
		lambda[term.Uidx].Add(p_exponent);
	}

	auto end = std::chrono::steady_clock::now();
	std::cout << "\ncounted " << counter << " tuples.\n";
	std::cout << std::chrono::duration <double, std::milli> (end-start).count() << " ms\n\n";

	return RET_OK;
}