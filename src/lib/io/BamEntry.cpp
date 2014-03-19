#include "BamEntry.hpp"

bool BamEntry::resolveCigar(uint32_t pos, PileupData& rv) const {
    int readPos = _cigarParser.getSnvPileupOffset(pos - start());
    if (readPos >= 0) {
        rv.tid = tid();
        rv.base = base(readPos);
        rv.quality = baseQuality(readPos);
        return true;
    }

    return false;
}
