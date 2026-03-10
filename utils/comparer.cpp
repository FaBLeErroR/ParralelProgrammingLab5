#include "comparer.h"
#include <limits>
#include <cmath>
#include <iostream>

bool compare(
    const std::vector<std::complex<double>> &original,
    const std::vector<std::complex<double>> &restored)
{
    const double EPS = std::numeric_limits<double>::epsilon();
    const double TOLERANCE = std::pow(2.0, -32);

    for (size_t i = 0; i < original.size(); i++)
    {
        double re_orig = original[i].real();
        double im_orig = original[i].imag();
        double re_rest = restored[i].real();
        double im_rest = restored[i].imag();

        double diff_re = std::abs(re_orig - re_rest);
        double diff_im = std::abs(im_orig - im_rest);

        if (re_orig == 0.0)
        {
            if (diff_re > TOLERANCE * EPS)
            {
                return false;
            }
        }
        else
        {
            double tolerance_re = TOLERANCE * std::abs(re_orig);
            if (diff_re > tolerance_re)
            {
                return false;
            }
        }

        if (im_orig == 0.0)
        {
            if (diff_im > TOLERANCE * EPS)
            {
                return false;
            }
        }
        else
        {
            double tolerance_im = TOLERANCE * std::abs(im_orig);
            if (diff_im > tolerance_im)
            {
                return false;
            }
        }
    }
    return true;
}