#pragma once

#include "io/BamEntry.hpp"
#include "Pileup.hpp"

#include <bam.h>
#include <deque>
#include <utility>
#include <cstdint>


class PileupBuffer {
public:
    PileupBuffer();
    ~PileupBuffer();

    uint32_t size() const;
    bool empty() const;

    const BamEntry* front() const {
        return _buf.front();
    }

    int tid() const;
    int32_t start() const;
    uint32_t end() const;

    Pileup* pileup(uint32_t pos) const;

    bool push(const BamEntry* b);
    void clear();
    void clearBefore(int tid, uint32_t pos);

    PosCompare cmp(const PileupBuffer& other, int32_t* xstart, uint32_t* xend) const;

protected:
    std::deque<const BamEntry*> _buf;
    uint32_t _end;
};

inline uint32_t PileupBuffer::size() const {
    return _buf.size();
}

inline bool PileupBuffer::empty() const {
    return _buf.empty();
}
