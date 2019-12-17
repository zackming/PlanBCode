#include <cmath>
#include "BlackScholes.h"


double blackScholesCall(const double& spot, const double& strike, const double& vol, const double& rate,
	const double& expiry, const double& div)
{
	auto d1 = (std::log(spot / strike) + (rate - div + 0.5 * vol * vol) * expiry) / (vol * std::sqrt(expiry));
	auto d2 = d1 - (vol * std::sqrt(expiry));
	auto callPrc = (spot * std::exp(-div * expiry) * normalCDF(d1)) - (strike * std::exp(-rate * expiry) * normalCDF(d2));

	return callPrc;
}

