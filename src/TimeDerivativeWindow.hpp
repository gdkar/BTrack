#pragma once

#include <fftw3.h>
#include "Math.hpp"
#include "Plan.hpp"

namespace retrack {
template<class InputIt, class OutputIt>
OutputIt time_derivative_window(InputIt wbegin, InputIt wend, OutputIt dst)
{
    auto n = int(std::distance(wbegin,wend));
    if(!n)
        return dst;
    auto time_r = std::vector<float>(wbegin,wend), time_i = std::vector<float>(n,0.f);
    {
        auto plan = Plan::dft_1d_c2c(n, &time_r[0],&time_i[0],&time_r[0],&time_i[0]);
        plan.execute();
        {
            auto base_idx = n/2 + 1;
            auto base_mul = - ( n - 1 ) * 0.5f;
            auto norm_mul = float(2 * M_PI) / (n * n);
            using std::swap;
            for(auto i = 0; i < n; ++i) {
                auto idx = (i + base_idx ) % n;
                auto mul =-(i + base_mul) * norm_mul;
                time_r[idx] *= mul;
                time_i[idx] *= mul;
                swap(time_r[idx],time_i[idx]);
            }
        }
        plan.execute();
    }
    return std::copy(time_r.cbegin(),time_r.cend(),dst);
}
}
