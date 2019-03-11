#pragma once

#include "sgstransversal.h"
#include "trivialtransversal.h"
#include "matrix.h"
#include "polynomial.h"

#include "utility.h"


using namespace CoffeeCode::Std;


template<typename IndexT>
struct MonomialAndIndex {
	using ExponentT = typename Polynomial::ExponentT;

	ExponentT exponent;
	IndexT Uidx;
};

template<typename MatrixT>
inline auto ChannelAction(const typename MatrixT::RowVectorT::StoreT XwoY, const typename MatrixT::RowVectorT::StoreT YwoX, const typename MatrixT::RowVectorT::StoreT XnY, const MatrixT& M)
{
	using VectorT = typename MatrixT::RowVectorT;

	const auto XwoYvec = VectorT(XwoY);
	const auto YwoXvec = VectorT(YwoX);
	const auto XnYvec = VectorT(XnY);

	// explicit cast allowed if ExponentT is large enough
	const auto u1 = XwoYvec.popcount();
	const auto u2 = YwoXvec.popcount();
	const auto u3 = XnYvec.popcount();

	const typename VectorT::StoreT Uidx = (M * (XwoYvec + YwoXvec) + (YwoXvec + XnYvec)).vec;

	return MonomialAndIndex<typename VectorT::StoreT>{
		Polynomial::MonomialT::MakeExponent(u1, u2, u3),
		Uidx
	};
}
template<typename MatrixT>
inline auto ChannelAction(const typename MatrixT::RowVectorT::StoreT subsetX, const typename MatrixT::RowVectorT::StoreT subsetY, const MatrixT& M)
{
	using StoreT = typename MatrixT::RowVectorT::StoreT;

	// this does integer promotion if StoreT is smaller than int,
	// so we implicitly cast back
	const auto XwoY = checked_cast<StoreT>(subsetX & ~subsetY);
	const auto YwoX = checked_cast<StoreT>(subsetY & ~subsetX);
	const auto XnY = checked_cast<StoreT>(subsetX & subsetY);

	return ChannelAction(XwoY, YwoX, XnY, M);
}


// print helper for accumulated lambdas
template<typename PolynomialT>
using ReducedLambdaT = std::unordered_map<
	PolynomialT,
	typename PolynomialT::CoefficientT,
	typename PolynomialT::Hash
>;
template<typename PolynomialT>
void PrintLambda(const ReducedLambdaT<PolynomialT>& lambda)
{
	size_t i = lambda.size();
	for (const auto& [poly, mult] : lambda) {
		std::cout << "[" << +mult << ", [" << poly << "]]";
		if (--i) std::cout << ",";
		std::cout << "\n";
	}
}