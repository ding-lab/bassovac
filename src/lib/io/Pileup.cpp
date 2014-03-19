#include "Pileup.hpp"
#include "utility/Lut.hpp"
#include "utility/CountingSort.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>

using namespace std;

Pileup* Pileup::create() {
    return new Pileup;
}

Pileup::Pileup() {
}

Pileup::~Pileup() {
}

void Pileup::set(std::vector<EntryType>& entries) {
    _entries.swap(entries);
}

uint32_t Pileup::readsMatching(int base, int minQual) const {
    uint32_t rv = 0;
    for (auto iter = _entries.begin(); iter != _entries.end(); ++iter) {
        if (iter->base == base && iter->quality >= minQual)
            ++rv;
    }
    return rv;
}

double Pileup::baseQualityHarmonicMean() const {
    double denominator = 0.0;
    for (auto iter = _entries.begin(); iter != _entries.end(); ++iter)
        denominator += Lut::phred2p_reciprocal(iter->quality);
    return size()/denominator;
}

vector<uint8_t> Pileup::baseQualities(int minQual) const {
    vector<uint8_t> qualities;
    qualities.reserve(_entries.size());
    for (auto iter = _entries.begin(); iter != _entries.end(); ++iter) {
        int qual = iter->quality;
        if (qual >= minQual) {
            qualities.push_back(qual);
            assert(qual < 256);
        }
    }
    countingSort<0, 255>(qualities.begin(), qualities.end());
    return qualities;
}

void Pileup::baseCounts(int bases[4]) const {
    memset(bases, 0, sizeof(int)*4);
    for (auto iter = _entries.begin(); iter != _entries.end(); ++iter) {
        if (iter->base & 1) ++bases[0];
        if (iter->base & 2) ++bases[1];
        if (iter->base & 4) ++bases[2];
        if (iter->base & 8) ++bases[3];
    }
}

int Pileup::variantAllele(int ref, int occ[4]) {
    int maxBase = -1;
    for (int i = 0; i < 4; ++i) {
        int value = 1 << i;
        if (ref & value || occ[i] < 1)
            continue;
        if (maxBase < 0 || occ[i] > occ[maxBase])
            maxBase = i;
    }

    if (maxBase < 0)
        return 0;
    return 1 << maxBase;
}
