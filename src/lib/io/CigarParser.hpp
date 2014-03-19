#pragma once

#include <cstdint>

class CigarParser {
public:
    struct PileupData {
        int tid;
        int base;
        uint8_t quality;
    };

    CigarParser(uint32_t const* cigar, int len);
    void advance();

    void advanceToNextRefBase();

    int getSnvPileupOffset(uint32_t pos);

    int currentOp() const {
        return currentOp_;
    }

    int currentOpLen() const {
        return currentOpLen_;
    }

    int currentOpIdx() const {
        return currentOpIdx_;
    }

    uint32_t refPos() const {
        return refPos_;
    }

    uint32_t readPos() const {
        return readPos_;
    }

private:
    uint32_t const* cigar_;
    int len_;
    uint32_t readPos_;
    uint32_t refPos_;
    int currentOpIdx_;
    int currentOp_;
    int currentOpLen_;
    bool started_;
};
