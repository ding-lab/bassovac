#pragma once

#include <bam.h>

#include <cstdint>


class BamFilter {
public:
    BamFilter()
        : _mask(0)
        , _minMapQual(0)
    {
    }

    BamFilter(uint32_t mask, uint32_t minMapQual)
        : _mask(mask)
        , _minMapQual(minMapQual)
    {
    }

    bool accept(const bam1_t* b) const {
        return (b->core.flag & _mask) == 0 && (b->core.qual >= _minMapQual);
    }

protected:
    uint32_t _mask;
    uint32_t _minMapQual;
};
