#include "fft_omp.h"
#include <vector>
#include <cmath>
#include <omp.h>

const double PI = std::acos(-1);

void fft_omp(std::complex<double>* output, std::size_t n,
             const std::complex<double>* input)
{
    if (n == 1)
    {
        output[0] = input[0];
        return;
    }

    std::size_t half = n / 2;

    std::vector<std::complex<double>> even(half);
    std::vector<std::complex<double>> odd(half);
    
    for (std::size_t i = 0; i < half; i++)
    {
        even[i] = input[2*i];
        odd[i]  = input[2*i+1];
    }

    std::vector<std::complex<double>> Fe(half);
    std::vector<std::complex<double>> Fo(half);

    #pragma omp task shared(Fe, even)
    fft_omp(Fe.data(), half, even.data());
    
    #pragma omp task shared(Fo, odd)
    fft_omp(Fo.data(), half, odd.data());
    
    #pragma omp taskwait

    for (std::size_t k = 0; k < half; k++)
    {
        std::complex<double> w = std::polar(1.0, -2.0 * PI * k / n);
        output[k] = Fe[k] + w * Fo[k];
        output[k + half] = Fe[k] - w * Fo[k];
    }
}

void ifft_omp(std::complex<double>* output, std::size_t n,
              const std::complex<double>* input)
{
    if (n == 1)
    {
        output[0] = input[0];
        return;
    }

    std::size_t half = n / 2;

    std::vector<std::complex<double>> even(half);
    std::vector<std::complex<double>> odd(half);
    
    for (std::size_t i = 0; i < half; i++)
    {
        even[i] = input[2*i];
        odd[i]  = input[2*i+1];
    }

    std::vector<std::complex<double>> Fe(half);
    std::vector<std::complex<double>> Fo(half);

    #pragma omp task shared(Fe, even)
    ifft_omp(Fe.data(), half, even.data());
    
    #pragma omp task shared(Fo, odd)
    ifft_omp(Fo.data(), half, odd.data());
    
    #pragma omp taskwait

    for (std::size_t k = 0; k < half; k++)
    {
        std::complex<double> w = std::polar(1.0, 2.0 * PI * k / n);
        output[k] = (Fe[k] + w * Fo[k]) / 2.0;
        output[k + half] = (Fe[k] - w * Fo[k]) / 2.0;
    }
}
