#pragma once

#include "types.h"


namespace CoffeeCode::NautyLink {
	using VColGCallbackType = std::function<void(const int*, size_t)>;

	void vcolg(const VColGCallbackType& callback, void *g, const size_t numcols, const size_t CallbackStride, const size_t CallbackOffset);
}