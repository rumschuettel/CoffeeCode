#include "CoffeeCode.h"

#include "ctmath.h"

#include <string>
#include <chrono>

using std::cout;
using std::endl;

#ifdef OPTIMIZE_FOR_DEPOLARIZING
using Polynomial = typename CoffeeCode::Polynomial<typename CoffeeCode::UnivariateMonomial>;
#else
using Polynomial = typename CoffeeCode::Polynomial<typename CoffeeCode::MultivariateMonomial<3>>;
#endif


namespace {
	enum ReturnValue {
		RET_OK = 0,
		RET_WRONG_INPUT = 1
	};

	// map reduce for lambdas
	// should match this format for use with CoffeeCode::PrintLambda
	[[maybe_unused]]
	auto ReduceLambda(const std::vector<Polynomial>& lambda)
	{
		ReducedLambdaT<Polynomial> out;

		// aggregate
		for (const auto& poly : lambda)
			out[poly] ++;
		
		return out;
	}
		
	// print helper for CCFull-specific output format
	[[maybe_unused]]
	void PrintLambda(const std::vector<Polynomial>& lambda)
	{
		size_t i = lambda.size();
		for (const auto& poly : lambda) {
			std::cout << "[1, [" << poly << "]]";
			if (--i) std::cout << ",";
			std::cout << "\n";
		}
	}
}

int FullSolver()
{
	using CoffeeCode::Bitmask;
	using CoffeeCode::BaseKSubsets;
	using CoffeeCode::ipow;

	constexpr auto K_TOT = K_SYS+K_ENV;

	// read matrix from cin
	std::string input = "";
	std::getline(std::cin, input);
	if (input.length() != (K_TOT * K_TOT)) {
		std::cerr << "wrong input size for K_SYS=" << K_SYS << " and K_ENV=" << K_ENV << endl;
		return RET_WRONG_INPUT;
	}

	const auto M = CoffeeCode::AdjacencyMatrix<K_SYS, K_ENV>::FromString(input);

	const auto start = std::chrono::steady_clock::now();

	// Calculate lambda and lambda_pre
	std::vector<Polynomial> lambda(ipow<size_t>(2, K_TOT)), lambda_pre(ipow<size_t>(2, K_SYS));

	// iterate over all U1, U2 and U3 subsets
	using VectorT = typename decltype(M)::RowVectorT;
	using SubsetT = typename VectorT::StoreT;

	// could overflow the type, hence we put the break condition at end
	for (SubsetT subsetX = 0; ; subsetX++) {
		for (SubsetT subsetY = 0; ; subsetY++) {
			const auto term = ChannelAction<Polynomial>(subsetX, subsetY, M);
			// extract system vertices; they're stored rightmost
			const SubsetT UAidx = term.Uidx & Bitmask<SubsetT, K_SYS>::mask0111;

			// add monomials
			lambda[term.Uidx] += term.exponent;
			lambda_pre[UAidx] += term.exponent;

			if (subsetY == BaseKSubsets<2, K_SYS>::count - 1) break;
		}

		if (subsetX == BaseKSubsets<2, K_SYS>::count - 1) break;
	}

	// Calculate lambda_a
	std::vector<Polynomial> lambda_a(ipow<size_t>(2, K_SYS));

	// build Ulookup
	const auto MAB = M.AB();
	using SubsetAT = typename decltype(MAB)::ColumnVectorT::StoreT;
	using SubsetBT = typename decltype(MAB)::RowVectorT::StoreT;
	for (SubsetAT subsetA = 0; ; subsetA++) {
		for (SubsetBT subsetB = 0; ; subsetB++) {
			const auto Ulookup = (MAB * subsetB + subsetA).vec;
			lambda_a[subsetA] += lambda_pre[Ulookup];

			if (subsetB == BaseKSubsets<2, K_ENV>::count - 1) break;
		}

		if (subsetA == BaseKSubsets<2, K_SYS>::count - 1) break;
	}

	// EXPORT AS JSON
	// some statistics
	std::cout << "{\n\"tuples\": " << BaseKSubsets<4, K_SYS>::count << ",\n\"time\": ";
	const auto end = std::chrono::steady_clock::now();
	std::cout << std::chrono::duration <double, std::milli> (end-start).count() << ",\n";

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