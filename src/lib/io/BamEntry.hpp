#pragma once

#include "CigarParser.hpp"

#include <bam.h>

#include <cstdint>

enum PosCompare {
    BEFORE,
    OVERLAP,
    AFTER
};

class BamEntry {
public:
    struct PileupData {
        int tid;
        int32_t readOffset;
        int base;
        int quality;

        bool operator==(PileupData const& rhs) const {
            return tid == rhs.tid
                && base == rhs.base
                && quality == rhs.quality
                ;
        }
    };

    explicit BamEntry(bam1_t* rawData);
    ~BamEntry();

    PosCompare cmp(const BamEntry& rhs) const;
    bool resolveCigar(uint32_t pos, PileupData& rv) const;

    const char* name() const;
    int32_t tid() const;
    int32_t start() const;
    // note: samtools uses uint32_t for end, int32_t for start
    uint32_t end() const;
    const bam1_t* rawData() const;
    int base(uint32_t offset) const;
    int baseQuality(uint32_t offset) const;
    uint32_t flag() const;
    uint32_t nCigar() const;

protected:
    bam1_t* _rawData;
    uint32_t _end;

    mutable uint32_t _readPos;
    mutable uint32_t _refPos;
    mutable int _lastCigarIdx;
    mutable CigarParser _cigarParser;
};

inline
BamEntry::BamEntry(bam1_t* rawData)
    : _rawData(rawData)
    , _end(bam_calend(&rawData->core, bam1_cigar(rawData)))
    , _readPos(0)
    , _refPos(0)
    , _lastCigarIdx(-1)
    , _cigarParser(bam1_cigar(rawData), rawData->core.n_cigar)
{
}

inline
BamEntry::~BamEntry() {
    bam_destroy1(_rawData);
}

inline
PosCompare BamEntry::cmp(const BamEntry& rhs) const {
    if (tid() < rhs.tid())
        return BEFORE;
    if (rhs.tid() < tid())
        return AFTER;

    if (end() <= uint32_t(rhs.start()))
        return BEFORE;
    if (rhs.end() <= uint32_t(start()))
        return AFTER;

    return OVERLAP;
}

inline
const char* BamEntry::name() const {
    return bam1_qname(_rawData);
}

inline
int32_t BamEntry::tid() const {
    return _rawData->core.tid;
}

inline
int32_t BamEntry::start() const {
    return _rawData->core.pos;
}

// note: samtools uses uint32_t for end, int32_t for start
inline
uint32_t BamEntry::end() const {
    return _end;
}

inline
const bam1_t* BamEntry::rawData() const {
    return _rawData;
}

inline
int BamEntry::base(uint32_t offset) const {
    return bam1_seqi(bam1_seq(_rawData), offset);
}

inline
int BamEntry::baseQuality(uint32_t offset) const {
    return bam1_qual(_rawData)[offset];
}

inline
uint32_t BamEntry::flag() const {
    return _rawData->core.flag;
}

inline
uint32_t BamEntry::nCigar() const {
    return _rawData->core.n_cigar;
}

