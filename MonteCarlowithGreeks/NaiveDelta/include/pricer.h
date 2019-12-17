#pragma once

#include <algorithm>
#include <array>
#include <functional>
#include <iostream>
#include <random>

double BlackScholesDelta(const double&, const double&, const double&,
	const double&, const double&, const double&, const double&);

double BlackScholesGamma(const double&, const double&, const double&,
	const double&, const double&, const double&, const double&);

std::vector<double> simulateNvPath(const double&, const double&, const double&,
	const double&, const double&, const size_t&);

std::tuple<double, double> simulateDeltaHedgedPath(const double&, const double&, const double&,
	const double&, const double&, const size_t&, const double&);

std::tuple<double, double, double> simulateDeltaGammaHedgedPath(const double&, const double&, const double&,
	const double&, const double&, const size_t&, const double&);

std::vector<double> simulateJumpPath(const double&, const double&, const double&,
	const double&, const double&, const size_t&, const double&, const double&,
	const double&, const double&);

std::tuple<double, int> simulateDeltaHedgedJumpPath(const double&, const double&, const double&,
	const double&, const double&, const size_t&, const double&,
	const double&, const double&, const double&, const double&);

std::tuple<double, double, double> simulateDeltaGammaHedgedJumpPath(const double&, const double&, const double&,
	const double&, const double&, const size_t&, const double&,
	const double&, const double&, const double&, const double&);




