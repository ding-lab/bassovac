#pragma once

#include "io/BamEntry.hpp"

#include <cstdint>
#include <vector>

class Pileup {
public:
    typedef BamEntry::PileupData EntryType;
    static int variantAllele(int ref, int baseCounts[4]);

    static Pileup* create();
    ~Pileup();

    void set(std::vector<EntryType>& entries);

    bool empty() const;
    uint32_t size() const;
    uint32_t readsMatching(int base, int minQual) const;
    double baseQualityHarmonicMean() const;

    const EntryType& operator[](uint32_t idx) const;

    std::vector<uint8_t> baseQualities(int minQual) const;

    void baseCounts(int occ[4]) const;

protected:
    Pileup();

protected:
    std::vector<EntryType> _entries;
};

inline bool Pileup::empty() const {
    return _entries.empty();
}

inline uint32_t Pileup::size() const {
    return _entries.size();
}

inline const Pileup::EntryType& Pileup::operator[](uint32_t idx) const {
    return _entries[idx];
}
