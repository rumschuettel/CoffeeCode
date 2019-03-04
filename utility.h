#pragma once

#include <bitset>

namespace {
	// debug helper
	template<typename T>
	void print(const T& vec)
	{
		for (auto it = vec.rbegin(); it != vec.rend(); it++)
			std::cout << +(*it) << " ";
		//std::cout << "\n";
	}

	template<size_t Len, typename T, typename BT = std::bitset<Len>>
	void print(const T& vec)
	{
		std::cout << BT(vec) << "\n";
	}
}