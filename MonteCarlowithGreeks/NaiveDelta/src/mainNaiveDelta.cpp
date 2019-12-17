#include "payoff.h"
#include "pricer.h"

int main)
{
	const size_t numSteps = 252;
	const size_t numReps = 100;
	auto spot = 41.0;
	auto strike = 40.0;
	auto rate = 0.08;
	auto vol = 0.30;
	auto div = 0.0;
	auto expiry = 1.0;
	auto callPrice = 0.0;
	auto dt = expiry / numSteps;
	auto errdt = std::exp((rate - div) * dt);
	auto beta1 = -1;
	auto CT = 0.0;
	auto sumCT = 0.0;
	auto sumCT2 = 0.0;

	for (auto i = 0; i < numReps; ++i)
	{
		auto result = simulateDeltaHedgedPath(spot, rate, vol, div, expiry, numSteps, strike);
		auto spotT = std::get<0>(result);
		auto cv = std::get<1>(result);
		callPrice += callPayoff(spotT, strike) + beta1 * cv;
		CT = callPayoff(spotT, strike) + beta1 * cv;
		sumCT += CT;
		sumCT2 += CT * CT;
	}

	callPrice = sumCT / numReps;
	callPrice *= std::exp(-rate * expiry);

	auto SD = std::sqrt((sumCT2 - (sumCT * sumCT) / numReps) * std::exp(-2 * rate * expiry) / (numReps - 1));
	auto SE = SD / sqrt(numReps);

	std::cout << "The Call price is: " << callPrice << std::endl;
	std::cout << "The Standard Error is: " << SE << std::endl;

	return 0;
}