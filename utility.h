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

template<typename B, typename A>
inline B checked_cast(const A a)
{
	B b = static_cast<B>(a);
	if (a != b) throw new std::exception();
	return b;
}

#else

#define checked_cast static_cast

#endif

