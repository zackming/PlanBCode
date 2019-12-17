#include "pricer.h"
#include "rand.h"
#include <numeric>
#include "cdfAndpdf.h"

std::vector<double> simulateNvPath(const double& spot, const double& rate, const double& vol,
	const double& div, const double& expiry, const size_t& numSteps)
{
	auto dt = expiry / numSteps;
	auto nudt = (rate - div - 0.5 * vol * vol) * dt;
	auto sigdt = vol * std::sqrt(dt);
	auto z = std::vector<double>(numSteps, 0.0);
	auto path = std::vector<double>(numSteps, spot);
	RandomNormals(z);

	for (auto i = 1; i < numSteps; ++i)
	{
		path[i] = path[i - 1] * std::exp(nudt + sigdt * z[i]);
	}

	return path;
}



double BlackScholesDelta(const double& spot, const double& t, const double& strike,
	const double& expiry, const double& vol, const double& rate, const double& div)
{
	auto tau = expiry - t;
	auto d1 = (std::log(spot / strike) + (rate - div + 0.5 * vol * vol) * tau) / (vol * std::sqrt(tau));
	auto delta = std::exp(-div * tau) * normalCDF(d1);

	return delta;
}

double BlackScholesGamma(const double& spot, const double& t, const double& strike,
	const double& expiry, const double& vol, const double& rate, const double& div)
{
	auto tau = expiry - t;
	auto d1 = (std::log(spot / strike) + (rate - div + 0.5 * vol * vol) * tau) / (vol * std::sqrt(tau));
	auto gamma = (std::exp(-div * tau) * normpdf(d1) / (spot * vol * std::sqrt(tau)));

	return gamma;
}

std::tuple<double, double> simulateDeltaHedgedPath(const double& spot, const double& rate, const double& vol,
	const double& div, const double& expiry, const size_t& numSteps, const double& strike)
{
	auto dt = expiry / numSteps;
	auto nudt = (rate - div - 0.5 * vol * vol) * dt;
	auto sigdt = vol * std::sqrt(dt);
	auto z = std::vector<double>(numSteps, 0.0);
	auto path = std::vector<double>(numSteps, spot);
	auto delta = std::vector<double>(numSteps, 0.0);
	RandomNormals(z);
	auto erddt = std::exp((rate - div) * dt);
	auto cv = 0.0;

	for (auto i = 1; i < numSteps; ++i)
	{
		path[i] = path[i - 1] * std::exp(nudt + sigdt * z[i]);
		auto t = i * dt;
		delta[i] = BlackScholesDelta(path[i - 1], t, strike, expiry, vol, rate, div);
		cv += delta[i] * (path[i] - path[i - 1] * erddt);
	}

	return std::make_tuple(path[numSteps - 1], cv);
}


std::tuple<double, double, double> simulateDeltaGammaHedgedPath(const double& spot, const double& rate, const double& vol,
	const double& div, const double& expiry, const size_t& numSteps, const double& strike)
{
	auto dt = expiry / numSteps;
	auto nudt = (rate - div - 0.5 * vol * vol) * dt;
	auto sigdt = vol * std::sqrt(dt);
	auto z = std::vector<double>(numSteps, 0.0);
	auto path = std::vector<double>(numSteps, spot);
	auto delta = std::vector<double>(numSteps, 0.0);
	auto gamma = std::vector<double>(numSteps, 0.0);
	RandomNormals(z);
	auto errdt = std::exp((rate - div) * dt);
	auto egamma = std::exp((2 * (rate - div) + vol * vol) * dt) - 2 * errdt + 1;
	auto cv1 = 0.0;
	auto cv2 = 0.0;

	for (auto i = 1; i < numSteps; ++i)
	{
		path[i] = path[i - 1] * std::exp(nudt + sigdt * z[i]);
		auto t = i * dt;
		delta[i] = BlackScholesDelta(path[i - 1], t, strike, expiry, vol, rate, div);
		gamma[i] = BlackScholesGamma(path[i - 1], t, strike, expiry, vol, rate, div);
		cv1 += delta[i] * (path[i] - path[i - 1] * errdt);
		cv2 += gamma[i] * ((path[i] - path[i - 1]) * (path[i] - path[i - 1]) - path[i - 1] * path[i - 1] * egamma);
	}

	return std::make_tuple(path[numSteps - 1], cv1, cv2);
}

std::vector<double> simulateJumpPath(const double& spot, const double& rate, const double& vol,
	const double& div, const double& expiry, const size_t& numSteps, const double& lamda, const double& k,
	const double& rate_j, const double& vol_j)
{
	auto dt = expiry / numSteps;
	auto nudt = (rate - div - 0.5 * vol * vol) * dt;
	auto sigdt = vol * std::sqrt(dt);
	auto z = std::vector<double>(numSteps, 0.0);
	auto path = std::vector<double>(numSteps, spot);
	auto m = std::vector<int>(numSteps, 0.0);
	auto sumw = 0.0;
	RandomNormals(z);
	RandomPoisson(m, 3 * dt);

	for (auto i = 1; i < numSteps; ++i)
	{
		auto w = std::vector<double>(m[i], 0.0);
		RandomNormals(w);
		auto sumw = std::accumulate(w.rbegin(), w.rend(), 0.0);
		auto theta1 = (rate - div - lamda * k - 0.5 * vol * vol) * dt + vol * std::sqrt(dt) * z[i];
		auto theta2 = m[i] * (rate_j - 0.5 * vol_j * vol_j) + vol_j * sumw;
		path[i] = path[i - 1] * std::exp(theta1) * std::exp(theta2);
	}

	return path;
}

std::tuple<double, int> simulateDeltaHedgedJumpPath(const double& spot, const double& rate, const double& vol,
	const double& div, const double& expiry, const size_t& numSteps, const double& strike,
	const double& lamda, const double& k, const double& rate_j, const double& vol_j)
{
	auto dt = expiry / numSteps;
	auto nudt = (rate - div - 0.5 * vol * vol) * dt;
	auto sigdt = vol * std::sqrt(dt);
	auto z = std::vector<double>(numSteps, 0.0);
	auto m = std::vector<int>(numSteps, 0.0);
	auto path = std::vector<double>(numSteps, spot);
	auto delta = std::vector<double>(numSteps, 0.0);
	RandomNormals(z);
	RandomPoisson(m, 3 * dt);
	auto errdt = std::exp((rate - div) * dt);
	auto cv = 0;

	for (auto i = 1; i < numSteps; ++i)
	{
		auto w = std::vector<double>(m[i], 0.0);
		RandomNormals(w);
		auto sumw = std::accumulate(w.rbegin(), w.rend(), 0.0);
		auto theta1 = (rate - div - lamda * k - 0.5 * vol * vol) * dt + vol * std::sqrt(dt) * z[i];
		auto theta2 = m[i] * (rate_j - 0.5 * vol_j * vol_j) + vol_j * sumw;
		path[i] = path[i - 1] * std::exp(theta1) * std::exp(theta2);
		auto t = i * dt;
		delta[i] = BlackScholesDelta(path[i - 1], t, strike, expiry, vol, rate, div);
		cv += delta[i] * (path[i] - path[i - 1] * errdt);
	}

	return std::make_tuple(path[numSteps - 1], cv);
}

std::tuple<double, double, double> simulateDeltaGammaHedgedJumpPath(const double& spot, const double& rate, const double& vol,
	const double& div, const double& expiry, const size_t& numSteps, const double& strike,
	const double& lamda, const double& k, const double& rate_j, const double& vol_j)
{
	auto dt = expiry / numSteps;
	auto nudt = (rate - div - 0.5 * vol * vol) * dt;
	auto sigdt = vol * std::sqrt(dt);
	auto z = std::vector<double>(numSteps, 0.0);
	auto m = std::vector<int>(numSteps, 0.0);
	auto path = std::vector<double>(numSteps, spot);
	auto delta = std::vector<double>(numSteps, 0.0);
	auto gamma = std::vector<double>(numSteps, 0.0);
	RandomNormals(z);
	RandomPoisson(m, 3 * dt);
	auto errdt = std::exp((rate - div) * dt);
	auto egamma = std::exp((2.0 * (rate - div) + vol * vol) * dt) - 2.0 * errdt + 1.0;
	auto cv1 = 0.0;
	auto cv2 = 0.0;

	for (auto i = 1; i < numSteps; ++i)
	{
		auto w = std::vector<double>(m[i], 0.0);
		RandomNormals(w);
		auto sumw = std::accumulate(w.rbegin(), w.rend(), 0.0);
		auto theta1 = (rate - div - lamda * k - 0.5 * vol * vol) * dt + vol * std::sqrt(dt) * z[i];
		auto theta2 = m[i] * (rate_j - 0.5 * vol_j * vol_j) + vol_j * sumw;
		path[i] = path[i - 1] * std::exp(theta1) * std::exp(theta2);
		auto t = i * dt;
		delta[i] = BlackScholesDelta(path[i - 1], t, strike, expiry, vol, rate, div);
		gamma[i] = BlackScholesGamma(path[i - 1], t, strike, expiry, vol, rate, div);
		cv1 += delta[i] * (path[i] - path[i - 1] * errdt);
		cv2 += gamma[i] * ((path[i] - path[i - 1]) * (path[i] - path[i - 1]) - path[i - 1] * path[i - 1] * egamma);
	}

	return std::make_tuple(path[numSteps - 1], cv1, cv2);
}