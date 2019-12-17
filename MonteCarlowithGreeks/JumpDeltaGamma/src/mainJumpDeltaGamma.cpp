#include "payoff.h"
#include "pricer.h"

int main()
{
	const size_t numSteps = 252;
	const size_t numReps = 100000;
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
	auto dt = expiry / numSteps;
	auto errdt = std::exp((rate - div) * dt);
	auto egamma = std::exp((2 * (rate - div) + vol * vol) * dt) - 2 * errdt + 1;
	auto beta1 = -1;
	auto beta2 = -0.5;
	auto CT = 0.0;
	auto sumCT = 0.0;
	auto sumCT2 = 0.0;

	for (auto i = 0; i < numReps; ++i)
	{
		auto result = simulateDeltaGammaHedgedJumpPath(spot, rate, vol, div, expiry, numSteps, strike, lamda, k, rate_j, vol_j);
		auto spotT = std::get<0>(result);
		auto cv1 = std::get<1>(result);
		auto cv2 = std::get<2>(result);
		CT = callPayoff(spotT, strike) + beta1 * cv1 + beta2 * cv2;
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