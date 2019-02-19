#pragma once

#include "sgstransversal.h"
#include "matrix.h"
#include "polynomial.h"


#ifndef K_SYS
#error "need to specify -DK_SYS=number of system qubits"
#endif

#ifndef K_ENV
#error "need to specify -DK_ENV=number of environment qubits"
#endif


template<typename IndexT>
struct MonomialAndIndex {
	CoffeeCode::Monomial::ExponentT u1, u2, u3;
	IndexT Uidx;
};

template<typename MatrixT>
inline auto ChannelAction(const typename MatrixT::RowVectorT::StoreT XwoY, const typename MatrixT::RowVectorT::StoreT YwoX, const typename MatrixT::RowVectorT::StoreT XnY, const MatrixT& M)
{
	using VectorT = typename MatrixT::RowVectorT;
	using ExponentT = CoffeeCode::Monomial::ExponentT;

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
	using StoreT = typename MatrixT::RowVectorT::StoreT;

	// this does integer promotion if StoreT is smaller than int,
	// so we implicitly cast back
	const auto XwoY = static_cast<StoreT>(subsetX & ~subsetY);
	const auto YwoX = static_cast<StoreT>(subsetY & ~subsetX);
	const auto XnY = static_cast<StoreT>(subsetX & subsetY);

	return ChannelAction(XwoY, YwoX, XnY, M);
}
