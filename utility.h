#pragma once

#include <bitset>


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


// checked integer conversion conversion
#ifdef DEBUG

#include <boost/numeric/conversion/cast.hpp>
#define checked_cast boost::numeric_cast

#else

#define checked_cast static_cast

#endif

