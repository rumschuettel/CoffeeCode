#pragma once

namespace CoffeeCode {
	// debug helper
	template<typename T>
	void print(const T& vec)
	{
		for (const auto e : vec)
			std::cout << (int)e << " ";
		std::cout << "\n";

	}
}