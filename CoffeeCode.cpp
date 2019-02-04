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

constexpr auto RET_OK = 0;
constexpr auto RET_WRONG_INPUT = 1;



void test() {
	using CoffeeCode::SGSTransversal;
	using CoffeeCode::SGSGenerator;
	using CoffeeCode::Group;
	using CoffeeCode::Permutation;

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

	using TupleT = sgs::TupleT<100>;
	TupleT vec{ 0, 0, 0, 1, 0, 0, 0 };
	print(vec);

	std::cout << std::boolalpha << sgs::IsCanonical<100>(vec) << "\n";

	size_t count = 0;
	for (const auto& tuple : sgs::TupleCosets<8>()) {
		if ((count++ % 10000) > 0)
			continue;
		print(tuple);
	}
}

int main()
{
	test();
	return 0;


	constexpr auto k_sys = K_SYS, k_env = K_ENV, k_tot = k_sys + k_env;

	// read matrix from cin
	std::string input = "";
	std::getline(std::cin, input);
	if (input.length() != (k_tot * k_tot)) {
		std::cerr << "wrong input size for k_sys=" << k_sys << " and k_env=" << k_env << endl;
		return RET_WRONG_INPUT;
	}

	CoffeeCode::AdjacencyMatrix<k_sys, k_env> M{ input };

	// Calculate lambda and lambda_pre
	constexpr auto max_exponent = k_tot;
	std::vector<Polynomial<max_exponent>> lambda(ipow(2, k_tot)), lambda_pre(ipow(2, k_sys));

	// iterate over all U1, U2 and U3 subsets
	// we do not cache as memory access is potentially slower than counting
	using VectorT = decltype(M)::RowVectorT;
	using SubsetT = VectorT::StoreT;

	// could reach numerical limits, hence break condition at end
	for (SubsetT subsetX = 0; ; subsetX++) {
		for (SubsetT subsetY = 0; ; subsetY++) {
			// we explicitly allow downcasting, so no {} as initializer; in fact the compiler complains if we do.
			// note that the types match, though, as subsetX is of type VectorT::StoreT, and the bit operations
			// keep it within the allowed range.
			const auto XwoY = VectorT(subsetX & ~subsetY);
			const auto YwoX = VectorT(subsetY & ~subsetX);
			const auto XnY = VectorT(subsetX & subsetY);

			const Monomial::ExponentT u1 = static_cast<Monomial::ExponentT>(XwoY.popcount());
			const Monomial::ExponentT u2 = static_cast<Monomial::ExponentT>(YwoX.popcount());
			const Monomial::ExponentT u3 = static_cast<Monomial::ExponentT>(XnY.popcount());

			const SubsetT Uidx = (M * (XwoY + YwoX) + (YwoX + XnY)).vec;
			const SubsetT UAidx = Uidx & Bitmask1s<k_sys>::mask;

			// add monomials
			assert(u1 + u2 + u3 < max_exponent);
			const Monomial::ExponentT p_exponent = u1 + u2 + u3;
  			lambda[Uidx].Add(p_exponent);
			lambda_pre[UAidx].Add(p_exponent);

			if (subsetY == BaseKSubsets<2, k_sys>::count-1) break;
		}

		if (subsetX == BaseKSubsets<2, k_sys>::count-1) break;
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

			if (subsetB == BaseKSubsets<2, k_env>::count-1) break;
		}

		if (subsetA == BaseKSubsets<2, k_sys>::count-1) break;
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
