#pragma once

#include "sgstransversal.h"
#include "matrix.h"
#include "polynomial.h"

template<typename IndexT>
struct MonomialAndIndex {
	Monomial<>::ExponentT u1, u2, u3;
	IndexT Uidx;
};

template<typename MatrixT>
inline auto ChannelAction(const typename MatrixT::RowVectorT::StoreT subsetX, const typename MatrixT::RowVectorT::StoreT subsetY, const MatrixT& M)
{
	using VectorT = typename MatrixT::RowVectorT;
	// cast down explicitly
	const auto XwoY = VectorT(subsetX & ~subsetY);
	const auto YwoX = VectorT(subsetY & ~subsetX);
	const auto XnY = VectorT(subsetX & subsetY);

	const Monomial<>::ExponentT u1 = static_cast<Monomial<>::ExponentT>(XwoY.popcount());
	const Monomial<>::ExponentT u2 = static_cast<Monomial<>::ExponentT>(YwoX.popcount());
	const Monomial<>::ExponentT u3 = static_cast<Monomial<>::ExponentT>(XnY.popcount());

	const typename VectorT::StoreT Uidx = (M * (XwoY + YwoX) + (YwoX + XnY)).vec;

	return MonomialAndIndex<typename VectorT::StoreT>{
		u1, u2, u3,
		Uidx
	};
}