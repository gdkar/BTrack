#pragma once

#include "Math.hpp"
#include "SysUtils.hpp"

namespace retrack {
template<typename T, typename S, typename W>
void cutShift(T *target,S *src,const W &window)
{
    auto targetSize = window.size();
    auto hws = targetSize / 2;
    auto wbeg = window.data();
    auto wmid = wbeg + hws;
    auto wend = wbeg + targetSize;
    auto smid = src + hws;
    auto tmid = target + (targetSize-hws);
    auto mult = [](auto x, auto y){return x*y;};
    std::transform(wmid,wend,smid,target,mult);
    std::transform(wbeg,wmid,src, tmid,  mult);
}
}

