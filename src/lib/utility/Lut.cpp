#include "Lut.hpp"

#include <cassert>
#include <cmath>
#include <memory>
#include <valarray>

using std::valarray;

namespace {
    valarray<double> _lgamma_lut;
    double _phred2p[256];
    double _phred2p_reciprocal[256];
    bool isInit = false;
    std::unique_ptr<Lut::RootsOfUnity<double>> _roots;

    void init_phred() {
        if (isInit) return;
        for (int i = 0; i < 256; ++i) {
            _phred2p[i] = pow(10, i/-10.0);
            _phred2p_reciprocal[i] = 1.0/_phred2p[i];
        }
        isInit = true;
    }

    void init_lgamma(uint32_t maxValue) {
        _lgamma_lut.resize(maxValue + 1);
        for (uint32_t i = 1; i <= maxValue; ++i)
            _lgamma_lut[i] = ::lgamma(i);
    }


}

namespace Lut {
    void init(uint32_t maxReadDepth /*= 50000*/) {
        init_phred();
        init_lgamma(maxReadDepth);
        _roots.reset(new RootsOfUnity<double>(maxReadDepth));
    }

    double lgamma(uint32_t x) {
        assert(x>0);
        if (x < _lgamma_lut.size())
            return _lgamma_lut[x];
        else
            return ::lgamma(x);
    }

    double const* lgamma_arr(uint32_t x) {
        assert(x>0 && x < _lgamma_lut.size());
        return &_lgamma_lut[x];
    }

    double phred2p(uint8_t phred) {
        assert(isInit);
        return _phred2p[phred];
    }

    double phred2p_reciprocal(uint8_t phred) {
        assert(isInit);
        return _phred2p_reciprocal[phred];
    }

    std::complex<double> const& rootsOfUnity(uint32_t n) {
        return (*_roots)(n);
    }

}
