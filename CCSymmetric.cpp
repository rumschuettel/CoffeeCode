#include "nautylink.h"

#include "CoffeeCode.h"
#include "utility.h"

#include <unordered_set>
#include <chrono>

#include <tuple>



namespace CoffeeCode {
	// extract compile time parameters
	template<typename T>
	struct SymmetricInstance {
		// parameters given
		constexpr static size_t K_TOT = K_SYS + K_ENV;

		// validate parameters
		static_assert(T::k_sys == K_SYS);
		static_assert(T::k_env == K_ENV);
		static_assert(T::adjacency_matrix.size() == K_TOT);
		static_assert(T::adjacency_matrix[0].size() == K_TOT);

		// adjacency matrix
		using MatrixT = CoffeeCode::AdjacencyMatrix<K_SYS, K_ENV>;
		constexpr static auto M = MatrixT{ T::adjacency_matrix };
		constexpr static auto MAB = M.AB();

		// types
		using RowVectorT = typename MatrixT::RowVectorT;
		using PolynomialT = typename T::Polynomial;

		// group functionality provider
	private:
		template<typename Q>
		using SymmetryProviderT = typename Q::sgs::template SymmetryProvider<MatrixT>;
		using SymmetryProvider = detected_or_t<
			NautyLink::NautyLink<MatrixT>,
			SymmetryProviderT,
			T
		>;
		
	public:
		constexpr static auto GroupLink()
		{
			return SymmetryProvider(M);
		}

		using CanonicalImageT = typename SymmetryProvider::CanonicalImageT;
		using CanonicalImageHashT = typename SymmetryProvider::CanonicalImageHashT;

		// lambdas as hashmaps
		using LambdaT = unordered_map<
			CanonicalImageT,
			std::pair<
				typename T::Polynomial,
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

	using CoffeeCode::AdjacencyMatrix;
	using CoffeeCode::TrivialSGSTransversal;
	using CoffeeCode::TrivialSGSOrbit;

	template<size_t S>
	using AdjacencyMatrixT = std::array<std::array<CoffeeCode::BitType, S>, S>;

	using UnivariatePolynomial = typename CoffeeCode::Polynomial<typename CoffeeCode::UnivariateMonomial>;
	using MultivariatePolynomial = typename CoffeeCode::Polynomial<typename CoffeeCode::MultivariateMonomial<3>>;
	template<auto& SamplePoints>
	using SampledPolynomial = typename CoffeeCode::Polynomial<typename CoffeeCode::SampledPolynomial<SamplePoints>>;
	template<size_t N>
	using UnivariateSamples = std::array<std::array<double, 1>, N>;
	template<size_t N>
	using MultivariateSamples = std::array<std::array<double, 3>, N>;
	#include "cc-instance.h"

	using instance = CoffeeCode::SymmetricInstance<graphstate_instance>;
	using PolynomialT = instance::PolynomialT;


	// for numerical polynomial types, we don't want to reduce the lambdas
	// since hashing is generally a bad idea for floats.
	enum {
		REDUCE_IDENTITY,
		REDUCE_SIMPLIFY,
		REDUCE_CALCULATE
	};
#ifdef REDUCE_LAMBDA_IF_POSSIBLE
	constexpr auto REDUCE_LAMBDA = (std::is_same_v<PolynomialT, UnivariatePolynomial> || std::is_same_v<PolynomialT, MultivariatePolynomial>) ?
		REDUCE_SIMPLIFY : REDUCE_CALCULATE;
#else
	constexpr auto REDUCE_LAMBDA = REDUCE_IDENTITY;
#endif

	// map reduce for lambdas
	// should match this format for use with CoffeeCode::PrintLambda
	template<typename LambdaT>
	auto ReduceLambda(const LambdaT& lambda)
	{
		// aggregate
		ReducedLambdaT<PolynomialT> out;
		for (const auto&[key, value] : lambda) {
			const auto&[poly, mult] = value;
			out[poly] += checked_cast<typename PolynomialT::CoefficientT>(mult);
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

	// calculate Shannon entropy for lambda
	// need the indirection because ::value might not exist
	template<typename LambdaT, typename Q>
	void PrintShannonEntropyProxy(const LambdaT& lambda, const double divisor)
	{
		decltype(Q::value) aggregate = { 0 };
		CoffeeCode::MultiplicityType<2> total_multiplicity = 0;

		// aggregate Shannon entropy
		for (const auto&[key, value] : lambda) {
			const auto&[poly, mult] = value;
			const auto& values = poly.coefficients.value;

			total_multiplicity += mult;
			for (size_t s = 0; s < values.size(); s ++)
				if (values[s] != 0.)
					aggregate[s] -= mult * values[s] * divisor * log2(values[s] * divisor);
		}

		std::cout << total_multiplicity << ", [";
		// print as comma-separated list
		size_t i = aggregate.size();
		auto oldprec = std::cout.precision();
		std::cout.precision(std::numeric_limits<decltype(Q::value)::value_type>::digits10 + 1);
		std::cout << std::scientific;
		for (const auto val : aggregate) {
			std::cout << std::setprecision(15) << val;
			if (--i) std::cout << ",";
		}
		std::cout << "]\n";
		std::cout << std::defaultfloat;
		std::cout.precision(oldprec);
	}
	template<typename LambdaT>
	void PrintShannonEntropy(const LambdaT& lambda, const double divisor = 1.)
	{
		PrintShannonEntropyProxy<LambdaT, PolynomialT::CoefficientArrayT>(lambda, divisor);
	}

	// extract tuple and multiplicity from iterator;
	// if not provided will call appropriate group functions automatically
	template<typename GroupT, typename IteratorT>
	auto TupleAndStabMult(GroupT& group, const IteratorT& it) noexcept
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

#ifdef PARALLELIZE
#include <omp.h>
#endif

#ifndef MEASURE_ALL
	#define MEASURE_FILTER1(sth)
#else
	#ifdef PARALLELIZE
	#error "cannot parallelize and measure timings within threads; no speedup"
	#endif
	#define MEASURE_FILTER1(sth) sth
#endif

int SymmetricSolver() {
	using VectorT = typename instance::RowVectorT;
	using SubsetT = typename VectorT::StoreT;
	using CoefficientT = typename PolynomialT::CoefficientT;

	// PERFORMANCE MEASURE
	auto now = std::chrono::high_resolution_clock::now;
	const auto time_total_start = now();
	const decltype(now()) time_temp;
	MEASURE_FILTER1((decltype(now() - now()) time_nauty_CCA, time_nauty_CCB, time_nauty_CCC, time_nauty_CCD;))
	size_t counter_channel = 0;
	size_t counter_ptrace = 0;

	// CHANNEL ACTION
	const auto time_channel_start = now();

	// fill lambdas
	instance::LambdaT lambda, lambda_pre;

	// graph group functions
	auto group = instance::GroupLink();

	// iterate over colorings of graph
#ifdef PARALLELIZE
#if (_MSC_VER || __INTEL_COMPILER)
#define OMP_SHAREDEFAULT shared
#else
#define OMP_SHAREDEFAULT none
#endif
	#pragma omp parallel default(OMP_SHAREDEFAULT) firstprivate(group) shared(lambda, lambda_pre, counter_channel, std::cout)
	{
	const size_t THREAD_COUNT = omp_get_num_threads();
	const size_t THREAD_ID = omp_get_thread_num();

	size_t counter_channel_ = 0;
	instance::LambdaT lambda_, lambda_pre_;

	// poor man's parallelization
	// since the bottleneck for this loop is generally not the coloring iterator, we simply have every thread loop
	// so that the i'th thread only handles the i'th, THREAD_COUNT + i'th, 2*THREAD_COUNT + i'th, ...
	#pragma omp for nowait
	for (int i = 0; i < checked_cast<int>(THREAD_COUNT); i++)
	group.Colorings<4>(THREAD_COUNT, THREAD_ID, [&](const auto& tuple, const auto orbitSize4, const size_t) -> void {
		counter_channel_++;
#else
	group.Colorings<4>(1, 0, [&](const auto& tuple, const auto orbitSize4, const size_t) -> void{
		counter_channel++;
#endif

		// calculate multiplicity of base 4 tuple

		// get low and high bit from tuples
		// note that SubsetT is such that the i'th index equals a bit shift i to the left
		SubsetT subsetX{ 0 }, subsetY{ 0 };
		for (size_t i = 0; i < K_SYS; i++) {
			const auto low_bit = tuple[i] & 0b01;
			const auto high_bit = (tuple[i] & 0b10) >> 1;

			CoffeeCode::OrBit(subsetX, !!low_bit, i);
			CoffeeCode::OrBit(subsetY, !!high_bit, i);
		}

		// apply channel
		const auto term = ChannelAction<PolynomialT>(subsetX, subsetY, instance::M);
		
		//// A: Add to lambda
		{
			// multiplicity of resulting base 2 tuple
			MEASURE_FILTER1(time_temp = now();)
			const auto[UIdxCanonical, orbitSize2] = group.CanonicalColoring(term.Uidx);
			MEASURE_FILTER1(time_nauty_CCA += now() - time_temp;)

			const CoefficientT coeff = checked_cast<CoefficientT>(orbitSize4 / orbitSize2);

			// accumulate polynomial with potentially pre-existing terms
			// note that mult is equal for U123 that map to the same UIdxCanonical
		#ifdef PARALLELIZE
			auto& [poly, mult] = lambda_[UIdxCanonical];
		#else
			auto& [poly, mult] = lambda[UIdxCanonical];
		#endif
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
		#ifdef PARALLELIZE
			auto& [poly, mult] = lambda_pre_[UIdxCanonical_pre];
		#else
			auto& [poly, mult] = lambda_pre[UIdxCanonical_pre];
		#endif
			poly.Add(term.exponent, coeff);
			mult = orbitSize2_pre;
		}
	});


#ifdef PARALLELIZE
	#pragma omp critical
	{
	// combine hashmaps into one
	for (const auto& [U, val] : lambda_) {
		const auto& [poly, mult] = val;
		auto& [acc_poly, acc_mult] = lambda[U];
		acc_poly += poly;
		acc_mult = mult;
	}
	for (const auto& [U, val] : lambda_pre_) {
		const auto& [poly, mult] = val;
		auto& [acc_poly, acc_mult] = lambda_pre[U];
		acc_poly += poly;
		acc_mult = mult;
	}

	// accumulate counter
	counter_channel += counter_channel_;

	} // #pragma omp critical
	} // #pragma omp parallel
#endif

	const auto time_channel = now() - time_channel_start;


	// PARTIAL TRACE
	using SubsetBT = decltype(instance::MAB)::RowVectorT::StoreT;
	using SubsetAT = decltype(instance::MAB)::ColumnVectorT::StoreT;
	instance::LambdaT lambda_a;

	const auto time_ptrace_start = now();
	group.Colorings<2>(1, 0, [&](const auto& tuple, const auto, const size_t) -> void {
		counter_ptrace++;

		// TupleT to SubsetT and HashT because that one can be different
		SubsetAT subsetA{ 0 };
		for (size_t i = 0; i < K_SYS; i++)
			CoffeeCode::OrBit(subsetA, !!tuple[i], i);

		MEASURE_FILTER1(time_temp = now();)
		const auto[Akey, __] = group.CanonicalColoring(subsetA);
		MEASURE_FILTER1(time_nauty_CCC += now() - time_temp;)

		MEASURE_FILTER1(time_temp = now();)

		// this inner loop is as in CCFull with the break condition at the end to avoid overflows
		for (SubsetBT subsetB = 0; ; subsetB++) {
			// new key, subset in A but cast to full subset on entire graph
			const auto BtoA = checked_cast<SubsetT>((instance::MAB * subsetB + subsetA).vec);

			// find canonical image for key
			const auto[key, keyMult] = group.CanonicalColoring(BtoA);

			//// C: Add to lambda_a
			const auto&[poly_pre, mult_pre] = lambda_pre[key];
			auto&[poly_a, mult_a] = lambda_a[Akey];

			assert(keyMult == mult_pre);  // must hold by definition
			poly_a += poly_pre;
			mult_a = keyMult;

			if (subsetB == CoffeeCode::BaseKSubsets<2, K_ENV>::count - 1) break;
		}

		MEASURE_FILTER1((time_nauty_CCD += (now() - time_temp) / CoffeeCode::BaseKSubsets<2, instance::K_ENV>::count;))
	});
	const auto time_ptrace = now() - time_ptrace_start;

	

	// EXPORT AS JSON
	// some statistics
	auto duration = [](auto t) { return std::chrono::duration<double, std::milli>(t).count(); };
	auto const print_time = [&](auto what, auto t) { std::cout << "\"" << what << "\": " << duration(t) << ",\n"; };
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
	if constexpr (REDUCE_LAMBDA == REDUCE_SIMPLIFY)
		PrintLambda(ReduceLambda(lambda));
	else if constexpr (REDUCE_LAMBDA == REDUCE_CALCULATE)
		PrintShannonEntropy(lambda);
	else
		PrintLambda(lambda);
	std::cout << "],\n";
	// lambda_a
	std::cout << "\"lambda_a\": [\n";
	if constexpr (REDUCE_LAMBDA == REDUCE_SIMPLIFY)
		PrintLambda(ReduceLambda(lambda_a));
	else if constexpr (REDUCE_LAMBDA == REDUCE_CALCULATE)
		PrintShannonEntropy(lambda_a, 1./pow(2., K_ENV));
	else
		PrintLambda(lambda_a);
	std::cout << "]\n}\n";

	return RET_OK;
}
