#pragma once

#include "SysUtils.hpp"

#include <set>

namespace retrack {

template<class InputIt, class OutputIt>
OutputIt time_weighted_window(InputIt wbegin, InputIt wend, OutputIt dst)
{
    auto n = int(std::distance(wbegin,wend));
    auto c = -0.5f * (n - 1);
    for(; wbegin != wend; (c += 1.0f))
        *dst++ = *wbegin++ * c;
    return dst;
}
}
