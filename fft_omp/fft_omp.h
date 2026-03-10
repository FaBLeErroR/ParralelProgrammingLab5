#pragma once
#include <complex>
#include <cstddef>

void fft_omp(std::complex<double>* output, std::size_t n,
             const std::complex<double>* input);

void ifft_omp(std::complex<double>* output, std::size_t n,
              const std::complex<double>* input);
