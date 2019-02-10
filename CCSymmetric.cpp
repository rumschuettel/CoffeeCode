#include "CoffeeCode.h"

#include <string_view>

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

	// this is a compile-time problem instance
	struct double_legged_graph {
		using sgs = SGSTransversal<
			SGSGenerator<2, Group<
			Permutation<0, 4, 5, 6, 1, 2, 3>
			>>,
			SGSGenerator<3, Group<
			Permutation<0, 1, 3, 2, 4, 5, 6>
			>>,
			SGSGenerator<6, Group<
			Permutation<0, 1, 2, 3, 4, 6, 5>
			>>
			>;

		constexpr static size_t k_sys = 7, k_env = 1;
		constexpr static auto adjacency_matrix = AdjacencyMatrix<k_sys, k_env>( 73, 176, 64, 64, 134, 8, 8, 128 );
	};

	// extract compile time parameters
	template<typename T>
	struct Instance {
		// parameters given
		constexpr static size_t k_sys = T::k_sys;
		constexpr static size_t k_env = T::k_env;
		constexpr static size_t k_tot = k_sys + k_env;
		using sgs = typename T::sgs;
		static_assert(T::sgs::Length == k_sys);

		// adjacency matrix
		constexpr static auto M = T::adjacency_matrix;
	};
}


int SymmetricSolver() {
	using instance = Instance<double_legged_graph>;

	for (const auto& tuple : instance::sgs::TupleCosets<4>()) {
	}

	return RET_OK;
}