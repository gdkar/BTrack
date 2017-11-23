_Pragma("once")

#include "Math.hpp"
#include "SysUtils.hpp"

namespace retrack {

class ReSpectrum {
public:
    using value_type = float;
    using vector_type = std::vector<value_type>;
    using size_type = vector_type::size_type;
    using difference_type = vector_type::difference_type;
    using pointer = vector_type::pointer;
    using const_pointer = vector_type::const_pointer;
    static constexpr int value_alignment = 32;
    static constexpr int item_alignment  = value_alignment / sizeof(value_type);

protected:
    int m_size{};
    int m_coef{m_size?(m_size / 2 + 1):0};
    int m_spacing{align_up<int>(m_coef, item_alignment)};
    int64_t m_when{0};
public:
    value_type  epsilon{};

    vector_type X               = vector_type( size_type(spacing()) * 2, value_type{});
    vector_type M               = vector_type( size_type(spacing())    , value_type{});
    vector_type Phi             = vector_type( size_type(spacing())    , value_type{});
    vector_type mag             = vector_type( size_type(spacing())    , value_type{});

    vector_type dM_dt           = vector_type( size_type(spacing())    , value_type{});
    vector_type dPhi_dt         = vector_type( size_type(spacing())    , value_type{});

    vector_type dM_dw           = vector_type( size_type(spacing())    , value_type{});
    vector_type dPhi_dw         = vector_type( size_type(spacing())    , value_type{});

    vector_type d2Phi_dtdw      = vector_type( size_type(spacing())    , value_type{});

    vector_type d2M_dtdw        = vector_type( size_type(spacing())    , value_type{});

    ReSpectrum(int _size = 0) : m_size(_size){}
    ReSpectrum(ReSpectrum && ) noexcept = default;
    ReSpectrum & operator = (ReSpectrum && ) noexcept = default;
    ReSpectrum(const ReSpectrum &) = default;
    ReSpectrum&operator=(const ReSpectrum &) = default;

    pointer X_real() { return &X[0];}
    const_pointer X_real() const { return &X[0];}

    pointer X_imag() { return &X[spacing()];}
    const_pointer X_imag() const { return &X[spacing()];}

    pointer mag_data() { return &mag[0];}
    const_pointer mag_data() const { return &mag[0];}

    pointer M_data() { return &M[0];}
    const_pointer M_data() const { return &M[0];}

    pointer Phi_data() { return &Phi[0];}
    const_pointer Phi_data() const { return &Phi[0];}

    pointer dM_dt_data() { return &dM_dt[0];}
    const_pointer dM_dt_data() const { return &dM_dt[0];}

    pointer dM_dw_data() { return &dM_dw[0];}
    const_pointer dM_dw_data() const { return &dM_dw[0];}

    pointer dPhi_dt_data() { return &dPhi_dt[0];}
    const_pointer dPhi_dt_data() const { return &dPhi_dt[0];}

    pointer dPhi_dw_data() { return &dPhi_dw[0];}
    const_pointer dPhi_dw_data() const { return &dPhi_dw[0];}

    pointer d2Phi_dtdw_data() { return &d2Phi_dtdw[0];}
    pointer d2M_dtdw_data() { return &d2M_dtdw[0];}

    const_pointer d2Phi_dtdw_data() const { return &d2Phi_dtdw[0];}
    const_pointer d2M_dtdw_data() const { return &d2M_dtdw[0];}

    int size() const
    {
        return m_size;
    }
    int coefficients() const
    {
        return m_coef;
    }
    int spacing() const
    {
        return m_spacing;
    }
    int64_t when() const
    {
        return m_when;
    }
    void set_when(int64_t _when)
    {
        m_when = _when;
    }
    void resize(int _size);
    void reset(int _size, int64_t _when)
    {
        resize(_size);
        set_when(_when);
    }
    void unwrapFrom(const ReSpectrum &o);
};
}
