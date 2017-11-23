#pragma once

#include <array>
#include <initializer_list>
#include <algorithm>
#include <numeric>
#include <utility>
#include <memory>
#include <cstdint>
#include <cmath>
#include <cstdlib>

#include "Math.hpp"
#include "Plan.hpp"
#include "SysUtils.hpp"
#include "ReFFT.hpp"
#include "Range.hpp"
#include "ReSpectrum.hpp"

namespace retrack {

struct ACF {
    using value_type = Plan::float_type;
    using float_type = Plan::float_type;
    using size_type  = Plan::size_type;
    using difference_type = Plan::difference_type;
    using iterator = value_type *;
    using const_iterator = const value_type *;
    using range_type = Range<iterator>;
    using const_range_type = Range<const_iterator>;
    constexpr ACF() = default;
    ACF(ACF &&) noexcept = default;
    ACF &operator = (ACF && ) noexcept = default;
    ACF(size_type n)
    : m_size{n}
    , m_fsize{m_size * 2}
    , m_coef{m_fsize ? (m_fsize/2u) + 1u : 0u}
    , m_spacing{m_coef ? align_up<size_type>(m_coef,16u) : 0}
    {
        if(m_fsize) {
            m_time  = detail::make_fftwf<value_type>(m_fsize);
            m_freq  = detail::make_fftwf<value_type>(m_spacing * 2);
            m_freq1 = detail::make_fftwf<value_type>(m_spacing * 2);
            m_r2c = Plan::dft_1d_r2c(m_fsize, m_time.get(), m_freq.get(),m_freq.get() + m_spacing);
            m_c2r = Plan::dft_1d_c2r(m_fsize, m_freq.get(), m_freq.get() + m_spacing, m_time.get());
        }
    }
    size_type size() const
    {
        return m_size;
    }
    size_type fsize() const
    {
        return m_fsize;
    }
    size_type coefficients() const
    {
        return m_coef;
    }
    size_type spacing() const
    {
        return m_spacing;
    }
    range_type time_data()
    {
        return { m_time.get(), m_time.get() + fsize()};
    }
    const_range_type time_data() const
    {
        return { m_time.get(), m_time.get() + fsize()};
    }
    std::array<range_type,2> freq_data()
    {
        return { range_type{ m_freq.get(), m_freq.get() + coefficients()},
                 range_type{ m_freq.get() + m_spacing, m_freq.get() + m_spacing + coefficients()} };
;
    }
   std::array<const_range_type,2> freq_data() const
    {
        return { const_range_type{ m_freq.get(), m_freq.get() + coefficients()},
                 const_range_type{ m_freq.get() + m_spacing, m_freq.get() + m_spacing + coefficients()}};
    }
    std::array<range_type,2> freq1_data()
    {
        return { range_type{ m_freq1.get(), m_freq1.get() + coefficients()},
                 range_type{ m_freq1.get() + m_spacing, m_freq1.get() + m_spacing + coefficients()}};
;
    }
   std::array<const_range_type,2> freq1_data() const
    {
        return { const_range_type{ m_freq1.get(), m_freq1.get() + coefficients()},
                 const_range_type{ m_freq1.get() + m_spacing, m_freq1.get() + m_spacing + coefficients()}};
    }
    template<class I, class O>
    O process(I src, O dst)
    {
        auto _t = time_data();
        auto _f = freq_data();
        auto _fr = _f[0].begin();
        auto _fi = _f[1].begin();
        std::fill(
            std::copy_n(
                src
            , m_size
            , _t.begin())
        , _t.end()
        , value_type{});

        m_r2c.execute(_t.first,_fr);

        std::transform(_f[0].begin(),_f[0].end(),_fi,_f[0].begin(),[](auto r,auto i){return r*r+i*i;});
        std::fill(_f[1].begin(),_f[1].end(),value_type{});

        m_c2r.execute(_fr,_t.first);

        auto base_mult = value_type{1}/m_fsize;
        auto count = value_type(m_size);
        for(auto it = _t.begin(), et = it + m_size; it != et;)
            *dst++ = *it++ * base_mult/(count--);

        return dst;
    }
    template<class I, class O>
    O process_lapped(I src, O dst)
    {
        auto _t = time_data();
        auto _f0 = freq_data();
        auto _f1 = freq_data();

        std::copy_n(src, fsize(), _t.begin());
        m_r2c.execute(_t.begin(),_f0[0].begin());

        std::fill_n(_t.begin(), size(), float_type{});
        m_r2c.execute(_t.begin(),_f1[0].begin());

        auto _fr0 = _f0[0].first;
        auto _fi0 = _f0[1].first;
        auto _fr1 = _f1[0].first;
        auto _fi1 = _f1[1].first;
        for(auto i = 0; i < coefficients(); ++i) {
            auto r0 = _fr0[i],i0=_fi0[i],r1=_fr1[i],i1=_fi1[i];
            _fr0[i] = r0*r1 + i0*i1;
            _fi0[i] = i0*r1 - r0*i1;
        }
//        std::transform(_f0.begin(),_f0.end(),_f1.begin(),_f0.begin(),[](auto x,auto y){return x * std::conj(y);});

        m_c2r.execute(_fr0,_t.begin());

        auto base_mult = value_type{1}/(fsize() * size());
        return std::transform(_t.begin(),_t.begin() + size(), dst,[=](auto x){return x * base_mult;});
    }
protected:
    size_type                       m_size{};
    size_type                       m_fsize{};
    size_type                       m_coef{};
    size_type                       m_spacing{};
    detail::fftwf_ptr<value_type>   m_time{};
    detail::fftwf_ptr<value_type> m_freq{};
    detail::fftwf_ptr<value_type> m_freq1{};
    Plan                            m_r2c{};
    Plan                            m_c2r{};
};
}
