#include "Math.hpp"
#include "SysUtils.hpp"
#include "ReSpectrum.hpp"

using namespace retrack ;
void ReSpectrum::resize(int _size)
{
    if(_size == size())
        return;
    m_size = _size;
    m_coef = m_size / 2 + 1;
    m_spacing = align_up<int>(m_coef, item_alignment);
    X.resize(m_spacing * 2);
    cexpr_for_each(
          [sz=size_type(m_spacing)](auto & item){
              if(item.size() != sz)
                item.resize(sz);
          }
        , mag
        , M
        , Phi
        , dM_dw
        , dPhi_dw
        , dM_dt
        , dPhi_dt
        , d2Phi_dtdw
        , d2M_dtdw
        );
}
void ReSpectrum::unwrapFrom(const ReSpectrum &o)
{
    if(!size() || (o.size() != size())) {
        return;
    }

    const auto _dp0 = dPhi_dt_data();
    const auto _dp1 = o.dPhi_dt_data();
    const auto _p0 = Phi_data();
    const auto _p1 = o.Phi_data();

    const auto _half_hop   = 0.5f*(when() - o.when());
    const auto _twice_unit = 4 * value_type(M_PI) / size();

    auto i = 0;
    for(; i < m_coef; ++i) {
        auto _unit = _twice_unit * i;
        auto _incr = (_unit - *(_dp0 + i) - *(_dp1 + i)) * _half_hop;
        auto _next = *(_p1 + i) + _incr;
        _next = std::isnan(_next) ? value_type{} : _next;
        auto _prev = *(_p0 + i);
        auto _gen  = _next + princarg(_prev-_next);
        _p0[i] = std::isnan(_gen) ? value_type{} : _gen;
    }
}
