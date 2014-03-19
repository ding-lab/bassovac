#pragma once

#include <sam.h>
#include <bam.h>

#include <string>

class BamEntry;
class BamFilter;

struct Region {
    int tid;
    int beg;
    int end;
};


class BamReaderBase {
public:
    BamReaderBase();
    virtual ~BamReaderBase();

    virtual bam_header_t* header() const = 0;
    virtual std::string const& path() const = 0;

    virtual Region const* region() const {
        return 0;
    }

    // items returned by take must be deleted by the caller
    BamEntry* take();

    // items returned by peek should not be deleted by the caller
    BamEntry* peek();

    char const* targetName(int tid) const;

    void setFilter(BamFilter* filter) {
        filter_ = filter;
    }

protected:
    virtual bool takeImpl(bam1_t* entry) = 0;

private:
    BamEntry* buf_;
    BamFilter* filter_;
};
