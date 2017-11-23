#pragma once

#include "Plan.hpp"
#include "SysUtils.hpp"
//#include "VectorOpsComplex.hpp"
#include "TimeWeightedWindow.hpp"
#include "TimeDerivativeWindow.hpp"
#include "TimeAlias.hpp"
#include "ReSpectrum.hpp"

namespace retrack {

struct ReFFT {
    using value_type = float;
    using vector_type = std::vector<value_type>;
    using size_type = vector_type::size_type;
    using difference_type = vector_type::difference_type;
    using pointer = vector_type::pointer;
    using const_pointer = vector_type::const_pointer;
    static constexpr int value_alignment = 32;
    static constexpr int item_alignment  = value_alignment / sizeof(value_type);

    int m_size{}; /// <- transform size
    int m_coef{m_size ? (m_size / 2 + 1) : 0};
    int m_spacing{align_up(m_coef, item_alignment)};

    value_type m_epsilon = std::pow(10.0f, -80.0f/20.0f);//std::numeric_limits<float>::epsilon();

    value_type          m_time_width{};
    value_type          m_freq_width{};
    vector_type m_h     = vector_type(size());  /// <- transform window
    vector_type m_Dh    = vector_type(size());  /// <- time derivative of window
    vector_type m_Th    = vector_type(size());  /// <- time multiplied window
    vector_type m_TDh   = vector_type(size());  /// <- time multiplied time derivative window;

    detail::fftwf_ptr<value_type> m_flat  = detail::make_fftwf<value_type>(size()); /// <- input windowing scratch space.
    detail::fftwf_ptr<value_type> m_split = detail::make_fftwf<value_type>(spacing()*2); /// <- frequency domain scratch space.

    detail::fftwf_ptr<value_type> m_X     = detail::make_fftwf<value_type>(spacing()*2); /// <- transform of windowed signal
    detail::fftwf_ptr<value_type> m_X_Dh  = detail::make_fftwf<value_type>(spacing()*2); /// <- transform of time derivative windowed signal
    detail::fftwf_ptr<value_type> m_X_Th  = detail::make_fftwf<value_type>(spacing()*2); /// <- transform of time multiplied window
    detail::fftwf_ptr<value_type> m_X_TDh = detail::make_fftwf<value_type>(spacing()*2); /// <- transform of time multiplied time derivative window.

    Plan             m_plan_r2c{};
    Plan             m_plan_c2r{};
    void _finish_set_window();
    void _finish_process(float *src, ReSpectrum & dst, int64_t _when);

    ReFFT() = default;

    ReFFT ( ReFFT && ) noexcept = default;
    ReFFT &operator = ( ReFFT && ) noexcept = default;

    void initPlans();
    ReFFT ( int _size)
    : m_size{_size} { if(_size) initPlans(); }

    template<class It>
    It setWindow(It wbegin, It wend)
    {
        auto wn = std::distance(wbegin,wend);
        if(wn < m_size) {
            std::fill(m_h.begin(),m_h.begin() + (m_size-wn)/2,0.f);
            std::fill(std::copy(wbegin,wend, m_h.begin() + (m_size-wn)/2),m_h.end(),0.f);
        }else{
            std::copy_n(wbegin,size(), m_h.begin());
            wbegin += wn;
        }
        _finish_set_window();
        return wbegin;
    }
    template<class It>
    ReFFT( It wbegin, It wend)
    : ReFFT(int(std::distance(wbegin,wend)))
    {
        setWindow(wbegin,wend);
    }
    static ReFFT Kaiser(int _size, float alpha);
    virtual ~ReFFT();

    template<class It>
    void process( It src, It send, ReSpectrum & dst, int64_t when = 0)
    {
        auto tsrc = static_cast<float*>(alloca(
            m_size * sizeof(float))),
             tsend = tsrc + m_size;
        std::fill(std::copy(src,send,tsrc),tsend,0.0f);
        _finish_process(tsrc,dst,when);
    }
    template<class It>
    void process( It src, ReSpectrum & dst, int64_t when = 0)
    {
        process(src,std::next(src,m_size),dst,when);
    }
    template<class It, class iIt>
    void inverse( It dst, iIt _M, iIt _Phi)
    {
        auto _spacing = spacing();
        auto _coef    = m_coef;
        auto norm = 0.5f /(float(size()));
        std::transform(
            _M
           ,_M + _coef
           , &m_X[0]
           , [norm](auto x){
                return norm * std::exp(x);
            });
        std::copy(
            _Phi
           ,_Phi+ _coef
           , &m_X[0] + _spacing
            );
        for(auto i = 0; i < _coef; ++i) {
            auto m = m_X[i];
            auto t = m_X[i+_spacing];
            auto s = std::sin(t);
            auto c = std::cos(t);
            m_split[i] = c * m;
            m_split[i + _spacing] = s * m;
        }
        m_plan_c2r.execute(&m_split[0], &m_flat[0]);
        std::rotate_copy(
            &m_flat[0]
           ,&m_flat[0] + (m_size/2)
           ,&m_flat[m_size]
           ,dst
            );
    }
    template<class I, class O>
    void inverseCepstral(O dst, I src)
    {
        std::transform(src, src + m_coef, &m_split[0], std::log<value_type>);
        std::fill_n(&m_split[spacing()], m_coef, 0.0f);
        m_plan_c2r.execute(&m_split[0], &m_flat[0]);
        std::copy(&m_flat[0],&m_flat[m_size], dst);
    }
    int spacing() const;
    int size() const;
    int coefficients() const;
    const_pointer   h_data() const;
    const_pointer   Dh_data() const;
    const_pointer   Th_data() const;
    const_pointer   TDh_data() const;
    value_type      time_width() const;
    value_type      freq_width() const;
};
}
