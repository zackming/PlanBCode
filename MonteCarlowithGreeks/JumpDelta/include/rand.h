#pragma once
#include <random>



//extern "C" void RandomNormals(double* data, size_t size);
void RandomNormals(std::vector<double>&);


//extern "C" void RandomPoisson(double lambda, int* data, size_t size);
void RandomPoisson(std::vector<int>&, double lambda);