#pragma once

// iterating over a graph's tuples

namespace CoffeeCode::NautyLink {
	// nauty header
	#include <nauty.h>
	
	template<size_t Colors>
	struct VColGImpl;

	template<size_t Colors>
	class NautyOrbitIterator {
		friend struct VColGImpl<Colors>;

		const graph* G;

	public:
		NautyOrbitIterator(const graph* G)
			: G(G)
		{
			ColorGraph();
		}

	private:
		void ColorGraph();
		void Callback(FILE*, graph*, int*, int, int)
		{
		}

	};
}