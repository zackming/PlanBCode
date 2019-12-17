#include <cmath>
#include "cdfAndpdf.h"
double normalCDF(double x)
{
	return std::erfc(-x / std::sqrt(2)) / 2;
}

const double ONE_OVER_SQRT_2PI = 0.39894228040143267793994605993438;

double normpdf(double x)
{
	return ONE_OVER_SQRT_2PI * std::exp(-0.5 * x * x);
}