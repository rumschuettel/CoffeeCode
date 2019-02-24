#include "CoffeeCode.h"

#include "ctmath.h"

#include <string>
#include <chrono>

using std::cout;
using std::endl;


namespace {
	enum ReturnValue {
		RET_OK = 0,
		RET_WRONG_INPUT = 1
	};

	constexpr auto k_sys = K_SYS, k_env = K_ENV, k_tot = k_sys + k_env;

		
	// print helpers
	template<typename LambdaT>
	void PrintLambda(const LambdaT& lambda)
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
	using CoffeeCode::Monomial;
	using CoffeeCode::Polynomial;
	using CoffeeCode::Bitmask;
	using CoffeeCode::BaseKSubsets;
	using CoffeeCode::ipow;

	// read matrix from cin
	std::string input = "";
	std::getline(std::cin, input);
	if (input.length() != (k_tot * k_tot)) {
		std::cerr << "wrong input size for k_sys=" << k_sys << " and k_env=" << k_env << endl;
		return RET_WRONG_INPUT;
	}

	const auto M = CoffeeCode::AdjacencyMatrix<k_sys, k_env>::FromString(input);

	auto start = std::chrono::steady_clock::now();

	// Calculate lambda and lambda_pre
	std::vector<Polynomial> lambda(ipow(2, k_tot)), lambda_pre(ipow(2, k_sys));

	// iterate over all U1, U2 and U3 subsets
	using VectorT = decltype(M)::RowVectorT;
	using SubsetT = VectorT::StoreT;

	// could overflow the type, hence we put the break condition at end
	for (SubsetT subsetX = 0; ; subsetX++) {
		for (SubsetT subsetY = 0; ; subsetY++) {
			const auto term = ChannelAction(subsetX, subsetY, M);
			// extract system vertices; they're stored rightmost
			const SubsetT UAidx = term.Uidx & Bitmask<SubsetT, k_sys>::mask0111;

			// add monomials
			const Monomial::ExponentT p_exponent = term.uSum();
			lambda[term.Uidx].Add(p_exponent);
			lambda_pre[UAidx].Add(p_exponent);

			if (subsetY == BaseKSubsets<2, k_sys>::count - 1) break;
		}

		if (subsetX == BaseKSubsets<2, k_sys>::count - 1) break;
	}

	// Calculate lambda_a
	std::vector<Polynomial> lambda_a(ipow(2, k_sys));

	// build Ulookup
	const auto MAB = M.AB();
	using SubsetAT = decltype(MAB)::ColumnVectorT::StoreT;
	using SubsetBT = decltype(MAB)::RowVectorT::StoreT;
	for (SubsetAT subsetA = 0; ; subsetA++) {
		for (SubsetBT subsetB = 0; ; subsetB++) {
			const auto Ulookup = (MAB * subsetB + subsetA).vec;
			lambda_a[subsetA] += lambda_pre[Ulookup];

			if (subsetB == BaseKSubsets<2, k_env>::count - 1) break;
		}

		if (subsetA == BaseKSubsets<2, k_sys>::count - 1) break;
	}

	// EXPORT AS JSON
	// some statistics
	std::cout << "{\n\"tuples\": " << BaseKSubsets<4, k_sys>::count << ",\n\"time\": ";
	auto end = std::chrono::steady_clock::now();
	std::cout << std::chrono::duration <double, std::milli> (end-start).count() << ",\n";

	// lambda
	std::cout << "\"lambda\": [\n";
	PrintLambda(lambda);
	std::cout << "],\n";
	// lambda_a
	std::cout << "\"lambda_a\": [\n";
	PrintLambda(lambda_a);
	std::cout << "]\n}\n";

	return RET_OK;
}