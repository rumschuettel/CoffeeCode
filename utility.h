#pragma once

namespace {
	// debug helper
	template<typename T>
	void print(const T& vec)
	{
		for (const auto e : vec)
			std::cout << (int)e << " ";
		std::cout << "\n";
	}

	template<size_t Len, typename T, typename BT = std::bitset<Len>>
	void print(const T& vec)
	{
		std::cout << BT(vec) << "\n";
	}
}