#include "fft_async.h"
#include <vector>
#include <cmath>
#include <future>

const double PI = std::acos(-1);

void fft_async(std::complex<double>* output, std::size_t n,
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

    // Ограничение глубины: создаём потоки только если размер больше 256
    if (half > 256)  // для n=1024: half=512 - создаём, half=256 - уже нет
    {
        auto f1 = std::async(std::launch::async, fft_async, 
                             Fe.data(), half, even.data());
        auto f2 = std::async(std::launch::async, fft_async, 
                             Fo.data(), half, odd.data());
        
        f1.get();
        f2.get();
    }
    else
    {
        // Для маленьких подзадач - последовательные вызовы
        fft_async(Fe.data(), half, even.data());
        fft_async(Fo.data(), half, odd.data());
    }

    for (std::size_t k = 0; k < half; k++)
    {
        std::complex<double> w = std::polar(1.0, -2.0 * PI * k / n);
        output[k] = Fe[k] + w * Fo[k];
        output[k + half] = Fe[k] - w * Fo[k];
    }
}

void ifft_async(std::complex<double>* output, std::size_t n,
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

    if (half > 256)
    {
        auto f1 = std::async(std::launch::async, ifft_async, 
                             Fe.data(), half, even.data());
        auto f2 = std::async(std::launch::async, ifft_async, 
                             Fo.data(), half, odd.data());
        
        f1.get();
        f2.get();
    }
    else
    {
        ifft_async(Fe.data(), half, even.data());
        ifft_async(Fo.data(), half, odd.data());
    }

    for (std::size_t k = 0; k < half; k++)
    {
        std::complex<double> w = std::polar(1.0, 2.0 * PI * k / n);
        output[k] = (Fe[k] + w * Fo[k]) / 2.0;
        output[k + half] = (Fe[k] - w * Fo[k]) / 2.0;
    }
}
