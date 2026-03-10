#pragma once
#include <complex>
#include <cstddef>

void fft_seq(std::complex<double>* output, std::size_t n,
             const std::complex<double>* input);

void ifft_seq(std::complex<double>* output, std::size_t n,
              const std::complex<double>* input);