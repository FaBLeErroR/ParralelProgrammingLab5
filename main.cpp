#include <iostream>
#include <vector>
#include <complex>
#include <random>
#include <fstream>
#include <chrono>
#include <iomanip>

#include <omp.h>

#include "fft_seq/fft_seq.h"
#include "fft_omp/fft_omp.h"
#include "fft_async/fft_async.h"
#include "utils/comparer.h"

/// проверка корректности для всех трёх реализаций
void check_correctness(
    const std::vector<std::complex<double>> &input,
    std::vector<std::complex<double>> &fft_out,
    std::vector<std::complex<double>> &ifft_out,
    size_t n)
{
    fft_seq(fft_out.data(), n, input.data());
    ifft_seq(ifft_out.data(), n, fft_out.data());
    std::cout << (compare(input, ifft_out) ? "Seq correct" : "Seq incorrect") << "\n";

    fft_omp(fft_out.data(), n, input.data());
    ifft_omp(ifft_out.data(), n, fft_out.data());
    std::cout << (compare(input, ifft_out) ? "OMP correct" : "OMP incorrect") << "\n";

    fft_async(fft_out.data(), n, input.data());
    ifft_async(ifft_out.data(), n, fft_out.data());
    std::cout << (compare(input, ifft_out) ? "Async correct" : "Async incorrect") << "\n";
}

int main()
{
    const size_t n = 1024;
    const int repeat_times = 20;

    std::vector<std::complex<double>> input(n);
    std::vector<std::complex<double>> fft_out(n);
    std::vector<std::complex<double>> ifft_out(n);

    std::mt19937 gen(std::chrono::steady_clock::now().time_since_epoch().count());
    std::uniform_real_distribution<double> dist(-1.0, 1.0);

    for (auto &x : input)
        x = {dist(gen), dist(gen)};

    int max_threads = omp_get_max_threads();

    std::ofstream csv("results.csv");
    csv << std::fixed << std::setprecision(8);
    csv << "threads,seq_time_ms,seq_ifft_time_ms,omp_time_ms,omp_ifft_time_ms,async_time_ms,async_ifft_time_ms\n";

    for (int t = 1; t <= max_threads; t++)
    {
        omp_set_num_threads(t);

        double seq_time = 0, seq_ifft_time = 0;
        double omp_time = 0, omp_ifft_time = 0;
        double async_time = 0, async_ifft_time = 0;

        for (int r = 0; r < repeat_times; r++)
        {
            auto start = std::chrono::high_resolution_clock::now();
            fft_seq(fft_out.data(), n, input.data());
            auto end = std::chrono::high_resolution_clock::now();
            seq_time += std::chrono::duration<double, std::milli>(end - start).count();

            start = std::chrono::high_resolution_clock::now();
            ifft_seq(ifft_out.data(), n, fft_out.data());
            end = std::chrono::high_resolution_clock::now();
            seq_ifft_time += std::chrono::duration<double, std::milli>(end - start).count();
        }
        seq_time /= repeat_times;
        seq_ifft_time /= repeat_times;

        for (int r = 0; r < repeat_times; r++)
        {
            auto start = std::chrono::high_resolution_clock::now();
            fft_omp(fft_out.data(), n, input.data());
            auto end = std::chrono::high_resolution_clock::now();
            omp_time += std::chrono::duration<double, std::milli>(end - start).count();

            start = std::chrono::high_resolution_clock::now();
            ifft_omp(ifft_out.data(), n, fft_out.data());
            end = std::chrono::high_resolution_clock::now();
            omp_ifft_time += std::chrono::duration<double, std::milli>(end - start).count();
        }
        omp_time /= repeat_times;
        omp_ifft_time /= repeat_times;

        for (int r = 0; r < repeat_times; r++)
        {
            auto start = std::chrono::high_resolution_clock::now();
            fft_async(fft_out.data(), n, input.data());
            auto end = std::chrono::high_resolution_clock::now();
            async_time += std::chrono::duration<double, std::milli>(end - start).count();

            start = std::chrono::high_resolution_clock::now();
            ifft_async(ifft_out.data(), n, fft_out.data());
            end = std::chrono::high_resolution_clock::now();
            async_ifft_time += std::chrono::duration<double, std::milli>(end - start).count();
        }
        async_time /= repeat_times;
        async_ifft_time /= repeat_times;

        csv << t << ","
            << seq_time << "," << seq_ifft_time << ","
            << omp_time << "," << omp_ifft_time << ","
            << async_time << "," << async_ifft_time << "\n";
        ;
    }

    check_correctness(input, fft_out, ifft_out, n);

    csv.close();
    return 0;
}
