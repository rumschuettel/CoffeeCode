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


	// map reduce for lambdas
	// should match this format for use with CoffeeCode::PrintLambda
	template<typename LambdaT>
	auto ReduceLambda(const LambdaT& lambda)
	{
		// ugly, yeah...
		using PolynomialT = typename LambdaT::value_type::second_type::first_type;
		ReducedLambdaT<PolynomialT> out;

		// aggregate
		for (const auto& [key, value] : lambda) {
			const auto& [poly, mult] = value;
			assert(mult < PolynomialT::MaxCoefficient);
			out[poly] += static_cast<typename PolynomialT::CoefficientT>(mult);
		}
		
		return out;
	}
	
	// print helpers
	template<typename LambdaT>
	void PrintLambda(const LambdaT& lambda)
	{
		size_t i = lambda.size();
		for (const auto& [key, value] : lambda) {
			const auto& [poly, mult] = value;
			std::cout << "[" << +mult << ", [" << poly << "]]";
			if (--i) std::cout << ",";
			std::cout << "\n";
		}
	}
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

	// PERFORMANCE MEASURE
	auto now = std::chrono::high_resolution_clock::now;
	auto time_total_start = now();
	decltype(now()) time_temp;
	decltype(now() - now()) time_nauty_GO, time_nauty_CCA, time_nauty_CCB, time_nauty_CCC, time_nauty_CCD;
	size_t counter_channel = 0;
	size_t counter_ptrace = 0;


	// CHANNEL ACTION

	// lambdas as hash maps
	// TODO: replace with https://github.com/greg7mdp/sparsepp
	using LambdaT = std::unordered_map<
		CanonicalImageT,
		std::pair<
			CoffeeCode::Polynomial,
			CoffeeCode::NautyLink::OrbitSizeT
		>,
		CanonicalImageT::Hash
	>;
	LambdaT lambda, lambda_pre;


	auto time_channel_start = now();
	for (const auto& tuple : instance::sgs::TupleCosets<4>()) {
		counter_channel++;

		// calculate multiplicity of base 4 tuple
		time_temp = now();
		nauty.SetColoring(tuple);
		const auto stabGroupOrder4 = nauty.GroupOrder();
		time_nauty_GO += now() - time_temp;
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
			time_temp = now();
			const auto[UIdxCanonical, stabGroupOrder2] = nauty.CanonicalColoring(term.Uidx);
			time_nauty_CCA += now() - time_temp;
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
			time_temp = now();
			const auto[UIdxCanonical_pre, stabGroupOrder2_pre] = nauty.CanonicalColoring(
				term.Uidx & CoffeeCode::Bitmask<decltype(term.Uidx), instance::k_sys>::mask0111
			);
			time_nauty_CCB += now() - time_temp;
			const auto orbitSize2_pre = fullGroupOrder / stabGroupOrder2_pre;
			const CoefficientT coeff = static_cast<CoefficientT>(orbitSize4 / orbitSize2_pre);

			// accumulate polynomial with same terms as above
			// multiplicity is 1
			auto& [poly, mult] = lambda_pre[UIdxCanonical_pre];
			poly.Add(p_exponent, coeff);
			mult = orbitSize2_pre;
		}
	}
	auto time_channel = now() - time_channel_start;


	// PARTIAL TRACE
	using SubsetBT = decltype(instance::MAB)::RowVectorT::StoreT;
	using SubsetAT = decltype(instance::MAB)::ColumnVectorT::StoreT;
	LambdaT lambda_a;

	auto time_ptrace_start = now();
	for (const auto& tuple : instance::sgs::TupleCosets<2>()) {
		counter_ptrace ++;

		// TupleT to SubsetT and HashT because that one can be different
		SubsetAT subsetA{ 0 };
		for (size_t i = 0; i < instance::k_sys; i++) 
			CoffeeCode::OrBit(subsetA, !!tuple[i], i);

		time_temp = now();
		const auto [Akey, _] = nauty.CanonicalColoring(subsetA);
		time_nauty_CCC += now() - time_temp;

		time_temp = now();

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
		
		time_nauty_CCD += (now() - time_temp)/CoffeeCode::BaseKSubsets<2, instance::k_env>::count;
	}
	auto time_ptrace = (now() - time_ptrace_start);


	// EXPORT AS JSON
	// some statistics
	auto duration = [](auto t) { return std::chrono::duration<double, std::milli>(t).count(); };
	auto print_time = [&](auto what, auto t) { std::cout << "\"" << what << "\": " << duration(t) << ",\n"; };
	std::cout << "{\n";
	std::cout << "\"channel\": " << counter_channel << ",\n";
	std::cout << "\"ptrace\": " << counter_ptrace << ",\n";
	print_time("total time [ms]", now() - time_total_start);
	print_time("channel time [ms]", time_channel);
	print_time("ptrace time [ms]", time_ptrace);
	print_time("nauty GroupOrder [ms]", time_nauty_GO);
	print_time("nauty CanonicalColoring A [ms]", time_nauty_CCA);
	print_time("nauty CanonicalColoring B [ms]", time_nauty_CCB);
	print_time("nauty CanonicalColoring C [ms]", time_nauty_CCC);
	print_time("nauty CanonicalColoring D [ms]", time_nauty_CCD);

	// lambda
	std::cout << "\"lambda\": [\n";
#ifdef REDUCE_LAMBDA
	PrintLambda(ReduceLambda(lambda));
#else
	PrintLambda(lambda);
#endif
	std::cout << "],\n";
	// lambda_a
	std::cout << "\"lambda_a\": [\n";
#ifdef REDUCE_LAMBDA
	PrintLambda(ReduceLambda(lambda_a));
#else
	PrintLambda(lambda_a);
#endif
	std::cout << "]\n}\n";

	return RET_OK;
}