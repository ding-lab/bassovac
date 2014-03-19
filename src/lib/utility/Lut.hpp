#pragma once

#include <cassert>
#include <cmath>
#include <complex>
#include <cstdint>
#include <vector>

namespace Lut {
    void init(uint32_t maxReadDepth = 5000);
    double phred2p(uint8_t phred);
    double phred2p_reciprocal(uint8_t phred);
    double lgamma(uint32_t x);
    double const* lgamma_arr(uint32_t x);
    std::complex<double> const& rootsOfUnity(uint32_t n);

    template<typename RealType>
    struct RootsOfUnity {
        typedef std::complex<RealType> ValueType;

        RootsOfUnity(unsigned maxVal)
            : maxVal(maxVal)
            , table(maxVal + 1, 1)
        {
            static const std::complex<RealType> ci(0, 1);
            for (unsigned i = 2; i <= maxVal; ++i) {
                auto base = ci * RealType(2.0 * M_PI / i);
                table[i] = std::exp(base);
            }
        }

        ValueType const& operator()(unsigned n) const {
            assert(n <= maxVal && n >= 0);
            return table[n];
        }

        unsigned maxVal;
        std::vector<ValueType> table;
    };
}
