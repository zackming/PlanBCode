#include "payoff.h"
#include "pricer.h"

int main()
{
	const size_t numSteps = 252;
	const size_t numReps = 10;
	auto spot = 41.0;
	auto strike = 40.0;
	auto rate = 0.08;
	auto vol = 0.30;
	auto div = 0.0;
	auto lamda = 3;
	auto rate_j = -0.02;
	auto vol_j = 0.05;
	auto k = std::exp(rate_j) - 1;
	auto expiry = 1.0;
	auto callPrice = 0.0;
	auto CT = 0.0;
	auto sumCT = 0.0;
	auto sumCT2 = 0.0;

	for (auto i = 0; i < numReps; ++i)
	{
		auto path = simulateJumpPath(spot, rate, vol, div, expiry, numSteps, lamda, k, rate_j, vol_j);
		callPrice += callPayoff(path[numSteps - 1], strike);
		sumCT = callPrice;
		CT = callPayoff(path[numSteps - 1], strike);
		sumCT2 += CT * CT;
	}

	callPrice /= numReps;
	callPrice *= std::exp(-rate * expiry);

	auto SD = std::sqrt((sumCT2 - sumCT * sumCT / numReps) * std::exp(-2 * rate * expiry) / (numReps - 1));
	auto SE = SD / sqrt(numReps);

	std::cout << "The Call price is: " << callPrice << std::endl;
	std::cout << "The Standard Error is: " << SE << std::endl;
	return 0;
}