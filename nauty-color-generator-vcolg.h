#pragma once

#include "types.h"


namespace CoffeeCode::NautyLink {
	using VColGCallbackType = std::function<void(const int*, size_t)>;

	void vcolg(const VColGCallbackType& callback, void *g, size_t numcols);
}