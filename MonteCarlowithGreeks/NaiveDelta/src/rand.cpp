#include "rand.h"
#include <algorithm>

extern "C" void RandomNormals(double* data, size_t size)
{
	thread_local static auto mtgen = std::mt19937{ std::random_device{}() };
	auto norm = std::normal_distribution<>{ 0.0, 1.0 };
	std::generate(data, data + size, [&norm]() { return norm(mtgen); });
}

void RandomNormals(std::vector<double>& data)
{
	RandomNormals(data.data(), data.size());
}

extern "C" void RandomPoisson(double lambda, int* data, size_t size)
{
	thread_local static auto mtgen = std::mt19937{ std::random_device{}() };
	auto norm = std::poisson_distribution<>{ lambda };
	std::generate(data, data + size, [&norm]() { return norm(mtgen); });
	
}

void RandomPoisson(std::vector<int>& data, double lambda)
{
	RandomPoisson(lambda, data.data(), data.size());
}