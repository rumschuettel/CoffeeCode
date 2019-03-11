#include "nautylink.h"

#include "CoffeeCode.h"
#include "utility.h"

#include <unordered_map>
#include <unordered_set>
#include <chrono>

#include <tuple>


namespace CoffeeCode {
	// extract compile time parameters
	template<typename T>
	struct SymmetricInstance {
		// parameters given
		constexpr static size_t K_TOT = K_SYS + K_ENV;
		using sgs = typename T::sgs;

		// validate parameters
		static_assert(T::sgs::Length == K_SYS);
		static_assert(T::adjacency_matrix.size() == K_TOT);
		static_assert(T::adjacency_matrix[0].size() == K_TOT);

		// adjacency matrix
		using MatrixT = CoffeeCode::AdjacencyMatrix<K_SYS, K_ENV>;
		constexpr static auto M = MatrixT{ T::adjacency_matrix };
		constexpr static auto MAB = M.AB();

		// types
		using RowVectorT = typename MatrixT::RowVectorT;

		// group functionality provider
	private:
		template<typename Q>
		using SymmetryProviderT = typename Q::template SymmetryProvider<MatrixT>;
		using SymmetryProvider = detected_or_t<
			NautyLink::NautyLink<MatrixT>,
			SymmetryProviderT,
			sgs
		>;
		
	public:
		constexpr static auto GroupLink()
		{
			return SymmetryProvider(M);
		}

		using CanonicalImageT = typename SymmetryProvider::CanonicalImageT;
		using CanonicalImageHashT = typename SymmetryProvider::CanonicalImageHashT;

		// lambdas as hashmaps
		// TODO: replace with https://github.com/greg7mdp/sparsepp
		using LambdaT = std::unordered_map<
			CanonicalImageT,
			std::pair<
				CoffeeCode::Std::Polynomial,
				MultiplicityType<2>
			>,
			CanonicalImageHashT
		>;
	};

}

namespace {
	enum ReturnValue {
		RET_OK = 0,
		RET_WRONG_INPUT = 1
	};

	using namespace CoffeeCode::Std;

	using CoffeeCode::SGSTransversal;
	using CoffeeCode::SGSGenerator;
	using CoffeeCode::Group;
	using CoffeeCode::Permutation;
	using CoffeeCode::AdjacencyMatrix;

	using CoffeeCode::TrivialSGSTransversal;
	using CoffeeCode::TrivialSGSOrbit;

	template<size_t S>
	using AdjacencyMatrixT = std::array<std::array<CoffeeCode::BitType, S>, S>;

	#include "cc-instance.h"

	using instance = CoffeeCode::SymmetricInstance<graphstate_instance>;


	// map reduce for lambdas
	// should match this format for use with CoffeeCode::PrintLambda
	template<typename LambdaT>
	auto ReduceLambda(const LambdaT& lambda)
	{
		ReducedLambdaT<Polynomial> out;

		// aggregate
		for (const auto& [key, value] : lambda) {
			const auto& [poly, mult] = value;
			out[poly] += checked_cast<typename Polynomial::CoefficientT>(mult);
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

	// extract tuple and multiplicity from iterator;
	// if not provided will call appropriate group functions automatically
	template<typename GroupT, typename IteratorT>
	auto TupleAndStabMult(GroupT& group, const IteratorT& it) 
	{
		if constexpr( is_pair<std::decay_t<IteratorT>> ) {
			// iterator returns tuple and multiplicity; just pass along
			return it;
		}
		else {
			// iterator just returns tuple. Calculate tuple's multiplicity
			group.SetColoring(it);
			return std::make_pair(it, group.template ColoringMultiplicity<4>());
		}
	}
}


// group library
#include "utility.h"

#ifndef MEASURE_ALL
#define MEASURE_FILTER1(sth)
#else
#define MEASURE_FILTER1(sth) sth
#endif

int SymmetricSolver() {
	using VectorT = typename instance::RowVectorT;
	using SubsetT = typename VectorT::StoreT;
	using CoefficientT = typename Polynomial::CoefficientT;

	auto group = instance::GroupLink();

	// PERFORMANCE MEASURE
	auto now = std::chrono::high_resolution_clock::now;
	auto time_total_start = now();
	decltype(now()) time_temp;
	MEASURE_FILTER1((decltype(now() - now()) time_nauty_CCA, time_nauty_CCB, time_nauty_CCC, time_nauty_CCD;))
	size_t counter_channel = 0;
	size_t counter_ptrace = 0;


	// CHANNEL ACTION
	instance::LambdaT lambda, lambda_pre;

	auto time_channel_start = now();
	for (const auto& it : instance::sgs::TupleCosets<4>()) {
		const auto& [tuple, orbitSize4] = TupleAndStabMult(group, it);

		counter_channel++;

		// calculate multiplicity of base 4 tuple

		// get low and high bit from tuples
		// note that SubsetT is such that the i'th index equals a bit shit i to the left
		SubsetT subsetX{ 0 }, subsetY{ 0 };
		for (size_t i = 0; i < K_SYS; i++) {
			const auto low_bit = tuple[i] & 0b01;
			const auto high_bit = (tuple[i] & 0b10) >> 1;

			CoffeeCode::OrBit(subsetX, !!low_bit, i);
			CoffeeCode::OrBit(subsetY, !!high_bit, i);
		}

		// apply channel
		const auto term = ChannelAction(subsetX, subsetY, instance::M);
		
		//// A: Add to lambda
		{
			// multiplicity of resulting base 2 tuple
			MEASURE_FILTER1(time_temp = now();)
			const auto[UIdxCanonical, orbitSize2] = group.CanonicalColoring(term.Uidx);
			MEASURE_FILTER1(time_nauty_CCA += now() - time_temp;)

			const CoefficientT coeff = checked_cast<CoefficientT>(orbitSize4 / orbitSize2);

			// accumulate polynomial with potentially pre-existing terms
			// note that mult is equal for U123 that map to the same UIdxCanonical
			auto& [poly, mult] = lambda[UIdxCanonical];
			poly.Add(term.exponent, coeff);
			mult = orbitSize2;
		}

		//// B: Add to lambda_pre
		{
			// we project out the environment for term.Uidx
			MEASURE_FILTER1(time_temp = now();)
			const auto[UIdxCanonical_pre, orbitSize2_pre] = group.CanonicalColoring(
				term.Uidx & CoffeeCode::Bitmask<decltype(term.Uidx), K_SYS>::mask0111
			);
			MEASURE_FILTER1(time_nauty_CCB += now() - time_temp;)

			const CoefficientT coeff = checked_cast<CoefficientT>(orbitSize4 / orbitSize2_pre);

			// accumulate polynomial with same terms as above
			// multiplicity is 1
			auto& [poly, mult] = lambda_pre[UIdxCanonical_pre];
			poly.Add(term.exponent, coeff);
			mult = orbitSize2_pre;
		}
	}
	auto time_channel = now() - time_channel_start;


	// PARTIAL TRACE
	using SubsetBT = decltype(instance::MAB)::RowVectorT::StoreT;
	using SubsetAT = decltype(instance::MAB)::ColumnVectorT::StoreT;
	instance::LambdaT lambda_a;

	auto time_ptrace_start = now();
	for (const auto& it : instance::sgs::TupleCosets<2>()) {
		const auto& [tuple, _] = TupleAndStabMult(group, it);
		counter_ptrace ++;

		// TupleT to SubsetT and HashT because that one can be different
		SubsetAT subsetA{ 0 };
		for (size_t i = 0; i < K_SYS; i++) 
			CoffeeCode::OrBit(subsetA, !!tuple[i], i);

		MEASURE_FILTER1(time_temp = now();)
		const auto [Akey, __] = group.CanonicalColoring(subsetA);
		MEASURE_FILTER1(time_nauty_CCC += now() - time_temp;)

		MEASURE_FILTER1(time_temp = now();)

		// this inner loop is as in CCFull with the break condition at the end to avoid overflows
		for (SubsetBT subsetB = 0; ; subsetB++) {
			// new key, subset in A but cast to full subset on entire graph
			const auto BtoA = checked_cast<SubsetT>((instance::MAB * subsetB + subsetA).vec);

			// find canonical image for key
			const auto [key, keyMult] = group.CanonicalColoring(BtoA);

			//// C: Add to lambda_a
			const auto& [poly_pre, mult_pre] = lambda_pre[key];
			auto& [poly_a, mult_a] = lambda_a[Akey];

			assert(keyMult == mult_pre);  // must hold by definition
			poly_a += poly_pre;
			mult_a = keyMult;
			
			if (subsetB == CoffeeCode::BaseKSubsets<2, K_ENV>::count - 1) break;
		}
		
		MEASURE_FILTER1((time_nauty_CCD += (now() - time_temp)/CoffeeCode::BaseKSubsets<2, instance::K_ENV>::count;))
	}
	auto time_ptrace = now() - time_ptrace_start;

	

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

	MEASURE_FILTER1(
		print_time("group CanonicalColoring A [ms]", time_nauty_CCA);
		print_time("group CanonicalColoring B [ms]", time_nauty_CCB);
		print_time("group CanonicalColoring C [ms]", time_nauty_CCC);
		print_time("group CanonicalColoring D [ms]", time_nauty_CCD);
	)


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