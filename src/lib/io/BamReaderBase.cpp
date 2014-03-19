#include "BamReaderBase.hpp"

#include <cassert>

#include "BamEntry.hpp"
#include "BamFilter.hpp"

BamReaderBase::BamReaderBase()
    : buf_(0)
    , filter_(0)
{
}

BamReaderBase::~BamReaderBase() {
    delete buf_;
    buf_ = 0;
}

BamEntry* BamReaderBase::peek() {
    if (!buf_)
        buf_ = take();
    return buf_;
}

BamEntry* BamReaderBase::take() {
    if (buf_) {
        BamEntry* rv = buf_;
        buf_ = NULL;
        return rv;
    }

    bam1_t* entry = bam_init1();
    memset(entry, 0, sizeof(bam1_t));
    while (takeImpl(entry)) {
        if (!filter_ || filter_->accept(entry)) {
            return new BamEntry(entry);
        }
    }
    bam_destroy1(entry);
    return 0;
}

char const* BamReaderBase::targetName(int tid) const {
    assert(tid >= 0 && tid < header()->n_targets);
    return header()->target_name[tid];
}
