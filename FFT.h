#pragma once
#include <complex>
#include <iostream>
#include "math.h"
#define _USE_MATH_DEFINE

const double pi = 3.14159265358979323846; 
const std::complex<double> i(0.0f, -1.0f);


std::complex<double>* dft(double* f, int N)
{
	std::complex<double>* coefficients = new std::complex<double>[N];
	if (N == 1)
	{
		coefficients[0] = std::complex<double>(f[0], 0.0);
		return coefficients;
	}
	
	for (int k = 0; k < N; k++)
	{
		coefficients[k] = (0.0,0.0);
		for (int n = 0; n < N; n++)
		{
			std::complex<double> C = (0.0, 0.0);
			C = f[n] * std::exp((-2 * pi * k * n * i) / (double)N);
			coefficients[k] += C;
		}
	}
	return coefficients;
}

std::complex<double>* fft(double* f, int N)
{
	std::complex<double>* coefficients = new std::complex<double>[N];
	if (N == 1)
	{
		coefficients[0] = std::complex<double>(f[0],0.0);
		return coefficients;
	}
	
	std::complex<double>* Even = new std::complex<double>[N/2];
	std::complex<double>* Odd = new std::complex<double>[N/2];

	double* f1 = new double[N / 2];
	double* f2 = new double[N / 2];

	for (int i = 0; i < N/2; i++)
	{
		Even[i] = f[2*i];
	}

	for (int i = 0; i <N/2; i++)
	{
		Odd[i] = f[2 * i + 1];
	}

	Even = fft(reinterpret_cast<double*>(Even),N/2);
	Odd = fft(reinterpret_cast<double*>(Odd), N/2);

	for (int k = 0; k < (N/2); k++)
	{
		coefficients[k] = Even[k] + exp((-2 * pi * k * i) / (double)N) * Odd[k];
		coefficients[k+N/2] = Even[k] - exp((-2 * pi * k * i) / (double)N) * Odd[k];
	}
	return coefficients;
}