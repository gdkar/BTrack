#pragma once

#include <fftw3.h>
#include <algorithm>
#include <memory>
#include <utility>
#include <functional>
#include <utility>
#include <cstdint>
#include <cstddef>

namespace retrack {

namespace detail {
struct fftwf_mem_dtor {
    void operator () ( void *p ) const { fftwf_free(p); }
};
template<class T>
using fftwf_ptr = std::unique_ptr<T[], fftwf_mem_dtor>;
template<class T>
fftwf_ptr<T> make_fftwf(std::size_t n) {
    return fftwf_ptr<T>(static_cast<T*>(fftwf_malloc(n * sizeof(T))));
}
};

struct Plan {
    using value_type        = fftwf_plan;
    using reference         = value_type&;
    using const_reference   = const value_type&;
    using difference_type   = std::ptrdiff_t;
    using size_type         = std::size_t;
    using float_type        = float;

    typedef void (r2c_exec)(const fftwf_plan, float*ri,float*ro,float*io);
    typedef void (c2c_exec)(const fftwf_plan, float*ri,float*ii,float*ro,float*io);
    value_type  m_d{0};
    union {
        r2c_exec *m_r2c{};
        c2c_exec *m_c2c;
    };
    difference_type m_off_in {0};
    difference_type m_off_out{0};


    operator value_type() const { return m_d; }
    value_type get() const      { return m_d; }
    bool operator !() const     { return !m_d; }
    operator bool() const       { return !!m_d; }
    value_type release() { return std::exchange(m_d, value_type{}); }
    void reset()
    {
        if(m_d)
            fftwf_destroy_plan(m_d);
        m_d = 0;
        m_off_in = 0;
        m_off_out = 0;
        m_r2c = nullptr;
    }
    constexpr Plan() noexcept = default;
    explicit constexpr Plan(fftwf_plan _d, r2c_exec* fn, difference_type off_in, difference_type off_out) noexcept
    : m_d{_d}
    , m_r2c{fn}
    , m_off_in{off_in}
    , m_off_out{off_out}
    {
    }
    explicit constexpr Plan(fftwf_plan _d, c2c_exec* fn, difference_type off_in, difference_type off_out   ) noexcept
    : m_d{_d}
    , m_c2c{fn}
    , m_off_in{off_in}
    ,m_off_out{off_out}
    {
    }
    constexpr Plan(Plan && o) noexcept
    : m_d{o.m_d}
    , m_r2c{o.m_r2c}
    , m_off_in{o.m_off_in}
    , m_off_out{o.m_off_out}
    {
        o.m_d = value_type{0};
        o.m_r2c = nullptr;
        o.m_off_in = o.m_off_out = 0;
    }
    Plan &operator=(Plan && o) noexcept
    {
        swap(*this, o);
        return *this;
    }
   ~Plan()
    {
        if(m_d) {
            fftwf_destroy_plan(m_d);
            m_d = value_type{0};
        }
        m_r2c = nullptr;
    }
    friend void swap(Plan &lhs, Plan &rhs) noexcept
    {
        using std::swap;
        swap(lhs.m_d,       rhs.m_d);
        swap(lhs.m_r2c,     rhs.m_r2c);
        swap(lhs.m_off_in,  rhs.m_off_in);
        swap(lhs.m_off_out, rhs.m_off_out);
    }
    static Plan dft_1d_r2c(int _n, float_type *_ti, float_type *_ro, float_type *_io, unsigned flags = FFTW_ESTIMATE)
    {
        auto dims = fftwf_iodim{ _n, 1, 1};
        auto off_out = _io - _ro;
        return Plan(fftwf_plan_guru_split_dft_r2c(
            1, &dims, 0, nullptr, _ti, _ro, _ro + off_out, flags)
            , &fftwf_execute_split_dft_r2c, 0, off_out);
    }
    static Plan dft_1d_c2r(int _n, float_type *_ri, float_type *_ii, float_type *_to, unsigned flags = FFTW_ESTIMATE)
    {
        auto dims = fftwf_iodim{ _n, 1, 1};
        auto off_in = _ii - _ri;
        return Plan(fftwf_plan_guru_split_dft_c2r(
            1, &dims, 0, nullptr, _ri, _ri + off_in, _to, flags)
            , &fftwf_execute_split_dft_c2r, off_in, 0);
    }
    static Plan dft_1d_c2c(int _n, float_type *_ri, float_type *_ii, float_type *_ro,float_type *_io, unsigned flags = FFTW_ESTIMATE)
    {
        auto dims = fftwf_iodim{ _n, 1, 1};
        auto off_in = _ii - _ri;
        auto off_out= _io - _ro;
        return Plan(fftwf_plan_guru_split_dft(
            1, &dims, 0, nullptr, _ri, _ri + off_in, _ro, _ro + off_out, flags)
            , &fftwf_execute_split_dft, off_in, off_out);
    }
    void execute(float_type *_in, float_type *_out)
    {
        if(!m_d)
            return;
//            throw std::runtime_error("cannot execute unplanned FFTs.");
        if(!_in || !_out)
            return;
//            throw std::runtime_error("you gonna segfault, dumbass!");
        if(m_off_out && !m_off_in) {
            if(m_r2c)
                (*m_r2c)(
                    m_d
                  , _in
                  , _out
                  , _out + m_off_out
                );
        }else if(m_off_in && !m_off_out) {
            if(m_r2c)
                (*m_r2c)(
                    m_d
                  , _in
                  , _in + m_off_in
                  , _out
                );
        }else{
            if(m_c2c)
                (*m_c2c)(
                    m_d
                  , _in
                  , _in + m_off_in
                  , _out
                  , _out + m_off_out
                );
        }
    }
    void execute()
    {
        if(m_d)
            fftwf_execute(m_d);
    }
};
}
