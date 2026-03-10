#pragma once
#include <complex>
#include <cstddef>

void fft_async(std::complex<double>* output, std::size_t n,
               const std::complex<double>* input);

void ifft_async(std::complex<double>* output, std::size_t n,
                const std::complex<double>* input);
