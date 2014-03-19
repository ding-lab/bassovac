#include "CigarParser.hpp"

#include <bam.h>

#include <cassert>

namespace {
    inline bool isRefBase(int op) {
        return op == BAM_CMATCH || op == BAM_CDEL;
    }

    // bam.h documents, but defines no constants for the result of
    // bam_cigar_type:
    // bit 1: consume query; bit 2: consume reference
    static int const BAM_CONSUME_QUERY = 1;
    static int const BAM_CONSUME_REFERENCE = 2;
}

CigarParser::CigarParser(uint32_t const* cigar, int len)
    : cigar_(cigar)
    , len_(len)
    , readPos_(0)
    , refPos_(0)
    , currentOpIdx_(0)
    , currentOp_(bam_cigar_op(*cigar))
    , currentOpLen_(bam_cigar_oplen(*cigar))
    , started_(false)
{
}

void CigarParser::advance() {
    int type = bam_cigar_type(currentOp_);
    if (type & BAM_CONSUME_REFERENCE) {
        refPos_ += currentOpLen_;
    }
    if (type & BAM_CONSUME_QUERY) {
        readPos_ += currentOpLen_;
    }

    ++currentOpIdx_;
    assert(currentOpIdx_ < len_);
    currentOp_ = bam_cigar_op(cigar_[currentOpIdx_]);
    currentOpLen_ = bam_cigar_oplen(cigar_[currentOpIdx_]);
}

void CigarParser::advanceToNextRefBase() {
    while (currentOpIdx_ < len_ && !isRefBase(currentOp_)) {
        advance();
    }
}

int CigarParser::getSnvPileupOffset(uint32_t pos) {
    if (!started_) {
        advanceToNextRefBase();
        started_ = true;
    }

    while (pos - refPos_ >= currentOpLen_) {
        advance();
        advanceToNextRefBase();
    }

    if (currentOpIdx_ >= len_ || currentOp_ != BAM_CMATCH) {
        return -1;
    }
    return readPos_ + (pos - refPos_);
}
