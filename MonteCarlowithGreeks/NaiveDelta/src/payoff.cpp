#include "payoff.h"

double callPayoff(const double& spot, const double& strike)
{
	return std::max(spot - strike, 0.0);
}

double putPayoff(const double& spot, const double& strike)
{
	return std::max(strike - spot, 0.0);
}

