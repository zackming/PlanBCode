
#include "BlackScholes.h"
#include <iostream>

int main()
{
	auto spot = 41.0;
	auto strike = 40.0;
	auto rate = 0.08;
	auto vol = 0.30;
	auto div = 0.0;
	auto expiry = 1.0;

	auto callPrc = blackScholesCall(spot, strike, vol, rate, expiry, div);

	std::cout << "The Call price is: " << callPrc << std::endl;
	return 0;
}

