#include "CoffeeCode.h"

#include <string>

using std::cout;
using std::endl;


#ifndef K_SYS
#error "need to specify -DK_SYS=number of system qubits"
#endif

#ifndef K_ENV
#error "need to specify -DK_ENV=number of environment qubits"
#endif

namespace {
	enum ReturnValue {
		RET_OK = 0,
		RET_WRONG_INPUT = 1
	};

	constexpr auto k_sys = K_SYS, k_env = K_ENV, k_tot = k_sys + k_env;
}


int FullSolver()
{
	// read matrix from cin
	std::string input = "";
	std::getline(std::cin, input);
	if (input.length() != (k_tot * k_tot)) {
		std::cerr << "wrong input size for k_sys=" << k_sys << " and k_env=" << k_env << endl;
		return RET_WRONG_INPUT;
	}

	const auto M = CoffeeCode::AdjacencyMatrix<k_sys, k_env>::FromString(input);

	// Calculate lambda and lambda_pre
	constexpr auto max_exponent = k_tot;
	std::vector<Polynomial<max_exponent>> lambda(ipow(2, k_tot)), lambda_pre(ipow(2, k_sys));

	// iterate over all U1, U2 and U3 subsets
	using VectorT = decltype(M)::RowVectorT;
	using SubsetT = VectorT::StoreT;

	// could overflow the type, hence we put the break condition at end
	for (SubsetT subsetX = 0; ; subsetX++) {
		for (SubsetT subsetY = 0; ; subsetY++) {
			const auto [ u1, u2, u3, Uidx ] = ChannelAction(subsetX, subsetY, M);
			// TODO this might be a bug. Don't we need to have UAidx = Uidx >> k_env & ...?
			const SubsetT UAidx = Uidx & Bitmask1s<k_sys>::mask;

			// add monomials
			assert(u1 + u2 + u3 < max_exponent);
			const Monomial<>::ExponentT p_exponent = u1 + u2 + u3;
			lambda[Uidx].Add(p_exponent);
			lambda_pre[UAidx].Add(p_exponent);

			if (subsetY == BaseKSubsets<2, k_sys>::count - 1) break;
		}

		if (subsetX == BaseKSubsets<2, k_sys>::count - 1) break;
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

			if (subsetB == BaseKSubsets<2, k_env>::count - 1) break;
		}

		if (subsetA == BaseKSubsets<2, k_sys>::count - 1) break;
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

	return RET_OK;
}