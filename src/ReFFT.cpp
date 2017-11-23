#include <thread>
#include <mutex>
#include "ReFFT.hpp"
#include "KaiserWindow.hpp"
using namespace retrack;
namespace detail {
struct _wisdom_reg {
    static std::once_flag _wisdom_once;
    _wisdom_reg(){
        std::call_once(_wisdom_once,[](){
            fftwf_init_threads();
            fftwf_make_planner_thread_safe();
            fftw_init_threads();
            fftw_make_planner_thread_safe();
        });
    }
   ~_wisdom_reg() = default;
};
/*static*/ std::once_flag _wisdom_reg::_wisdom_once{};
_wisdom_reg the_registrar{};
}
void ReFFT::initPlans()
{
    if(m_size) {
        m_coef = m_size ? (m_size / 2 + 1) : 0; /// <- number of non-redundant forier coefficients
        m_spacing = align_up(m_coef, item_alignment);
        cexpr_for_each([sz=size_type(m_size)](auto & item){item.resize(sz);}
            , m_h
            , m_Dh
            , m_Th
            , m_TDh
              );
        m_flat = detail::make_fftwf<value_type>(m_size);
        cexpr_for_each([sz=size_type(spacing() * 2)](auto & item){item = detail::make_fftwf<value_type>(sz);}
            , m_split
            , m_X
            , m_X_Dh
            , m_X_Th
            , m_X_TDh
              );
        auto _real = &m_split[0]; auto _imag = &m_split[m_spacing]; auto _time = &m_flat[0];
        m_plan_r2c = Plan::dft_1d_r2c(m_size, _time, _real, _imag);
        m_plan_c2r = Plan::dft_1d_c2r(m_size, _real, _imag, _time);
    }
    _finish_set_window();
}
/*static*/ ReFFT ReFFT::Kaiser(int _size, float alpha)
{
    auto win    = vector_type(_size, 0.0f);
    make_kaiser_window(win.begin(),win.end(), alpha);
    return ReFFT(win.cbegin(),win.cend());
}

ReFFT::~ReFFT() = default;

void ReFFT::_finish_set_window()
{
    if(!m_size){
        m_time_width = 0.0f;
        m_freq_width = 0.0f;
    }else{
        time_weighted_window(m_h.cbegin(),m_h.cend(),m_Th.begin());
        time_derivative_window(m_h.cbegin(),m_h.cend(),m_Dh.begin());
        time_weighted_window(m_Dh.cbegin(),m_Dh.cend(),m_TDh.begin());

        auto norm_ = std::inner_product(&m_h[0], &m_h[0] + m_size, &m_h[0],value_type{});
        auto var_t_unorm = std::inner_product(&m_Th[0], &m_Th[0] + m_size, &m_Th[0],value_type{});
        auto var_w_unorm = std::inner_product(&m_Dh[0], &m_Dh[0] + m_size, &m_Dh[0],value_type{});

        m_time_width = 2 * std::sqrt(value_type(M_PI) * var_t_unorm) * (1.f/std::sqrt(norm_));
        m_freq_width = 2 * std::sqrt(value_type(M_PI) * var_w_unorm) * (1.f/std::sqrt(norm_));
    }

}
void ReFFT::_finish_process(float *src, ReSpectrum & dst, int64_t _when )
{
    {
        auto do_window = [&](auto &w, auto &v) {
            cutShift(&m_flat[0], src, w);
            m_plan_r2c.execute(&m_flat[0], &v[0]);
        };
        do_window(m_h , m_X    );
        do_window(m_Dh, m_X_Dh );
        do_window(m_Th, m_X_Th );
        do_window(m_TDh,m_X_TDh);
    }

    using std::tie; using std::make_pair; using std::copy; using std::get;

    auto i = 0;

    const auto _real = &m_X[0], _imag    = &m_X[m_spacing]
        ,_real_Dh = &m_X_Dh[0], _imag_Dh = &m_X_Dh[m_spacing]
        ,_real_Th = &m_X_Th[0], _imag_Th = &m_X_Th[m_spacing]
        ,_real_TDh = &m_X_TDh[0], _imag_TDh = &m_X_TDh[m_spacing]
        ;
    auto _cmul = [](auto r0, auto i0, auto r1, auto i1) {
        return make_pair(r0 * r1 - i0 * i1, r0 * i1 + r1 * i0);
        };
    auto _pcmul = [&](auto x0, auto x1) {
        return _cmul(get<0>(x0),get<1>(x0),
                     get<0>(x1),get<1>(x1));
    };

    dst.reset(m_size, _when);
    for(; i < m_coef; ++i) {
        auto _X_r = *(_real + i), _X_i = *(_imag + i);
        *(dst.X_real() + i) = _X_r;
        *(dst.X_imag() + i) = _X_i;
        {
            auto _X_mag = std::hypot(_X_i,_X_r);
            auto _X_phi = std::atan2(_X_i,_X_r);
            *(dst.mag_data() + i) = _X_mag;
            *(dst.M_data() + i) = std::log(_X_mag);
            *(dst.Phi_data() + i) = std::isnan(_X_phi) ? value_type{} : _X_phi;
        }
    }
    auto max_mag = *std::max_element(dst.mag_data(),dst.mag_data() + m_coef);
    if(!max_mag)
        max_mag = 1.0f;
    dst.epsilon = max_mag * m_epsilon;
    auto _cinv = [e=sqr(dst.epsilon)](auto r, auto i) {
        auto n = r*r + i*i;
        auto m = (n<e) ? value_type{} : value_type{1}/n;
        return make_pair(r * m , -i * m);
    };
    i = 0;
    auto store = [](auto && x, auto y) { *y = x;};
    for(; i < m_coef; ++i) {
        auto _X_r = *(_real + i), _X_i = *(_imag + i);
        tie(_X_r, _X_i) = _cinv(_X_r,_X_i);

        auto _Dh_over_X = _cmul( *(_real_Dh + i),*(_imag_Dh + i) ,_X_r, _X_i );

        store(get<0>(_Dh_over_X), &dst.dM_dt  [0] + i);
        store(get<1>(_Dh_over_X), &dst.dPhi_dt[0] + i);

        auto _Th_over_X = _cmul( *(_real_Th + i),*(_imag_Th + i) ,_X_r, _X_i );

        store(-get<1>(_Th_over_X), &dst.dM_dw  [0] + i);
        store( get<0>(_Th_over_X), &dst.dPhi_dw[0] + i);

        auto _TDh_over_X    = _cmul( *(_real_TDh + i),*(_imag_TDh + i) ,_X_r, _X_i );
        auto _Th_Dh_over_X2 = _pcmul(_Th_over_X,_Dh_over_X);

        store(get<0>(_TDh_over_X )  - get<1>(_Th_Dh_over_X2),&dst.d2Phi_dtdw[0] + i);
        store(-get<1>(_TDh_over_X ) + get<1>(_Th_Dh_over_X2),&dst.d2M_dtdw[0] + i);
    }
}
int ReFFT::spacing() const
{
    return m_spacing;
}
ReFFT::const_pointer ReFFT::h_data() const { return &m_h[0];}
ReFFT::const_pointer ReFFT::Dh_data() const { return &m_Dh[0];}
ReFFT::const_pointer ReFFT::Th_data() const { return &m_Th[0];}
ReFFT::const_pointer ReFFT::TDh_data() const { return &m_TDh[0];}
int ReFFT::size() const
{
    return m_size;
}
int ReFFT::coefficients() const
{
    return m_coef;
}
ReFFT::value_type ReFFT::time_width() const
{
    return m_time_width;
}
ReFFT::value_type ReFFT::freq_width() const
{
    return m_freq_width;
}
