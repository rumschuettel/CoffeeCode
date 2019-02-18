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
inline auto ChannelAction(const typename MatrixT::RowVectorT::StoreT XwoY, const typename MatrixT::RowVectorT::StoreT YwoX, const typename MatrixT::RowVectorT::StoreT XnY, const MatrixT& M)
{
	using VectorT = typename MatrixT::RowVectorT;
	using ExponentT = Monomial<>::ExponentT;

	const auto XwoYvec = VectorT(XwoY);
	const auto YwoXvec = VectorT(YwoX);
	const auto XnYvec = VectorT(XnY);

	// explicit cast allowed if ExponentT is large enough
	const auto u1 = static_cast<ExponentT>(XwoYvec.popcount());
	const auto u2 = static_cast<ExponentT>(YwoXvec.popcount());
	const auto u3 = static_cast<ExponentT>(XnYvec.popcount());

	const typename VectorT::StoreT Uidx = (M * (XwoYvec + YwoXvec) + (YwoXvec + XnYvec)).vec;

	return MonomialAndIndex<typename VectorT::StoreT>{
		u1, u2, u3,
			Uidx
	};
}
template<typename MatrixT>
inline auto ChannelAction(const typename MatrixT::RowVectorT::StoreT subsetX, const typename MatrixT::RowVectorT::StoreT subsetY, const MatrixT& M)
{
	// cast down explicitly
	const auto XwoY = subsetX & ~subsetY;
	const auto YwoX = subsetY & ~subsetX;
	const auto XnY = subsetX & subsetY;

	return ChannelAction(XwoY, YwoX, XnY, M);
}
