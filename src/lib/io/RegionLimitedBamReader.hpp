#include "BamReaderBase.hpp"

#include <boost/format.hpp>
#include <stdexcept>
#include <string>

class RegionLimitedBamReader : public BamReader {
public:
    RegionLimitedBamReader(std::string const& path, char const* region);
    ~RegionLimitedBamReader();

    Region const* region() const;

protected:
    bool takeImpl(bam1_t* entry);

protected:
    std::string regionString_;
    bam_index_t* index_;
    bam_iter_t iter_;
    Region region_;
};

inline
RegionLimitedBamReader::RegionLimitedBamReader(std::string const& path, char const* region)
    : BamReader(path)
    , regionString_(region)
    , index_(bam_index_load(path.c_str()))
{
    using boost::format;
    if (!index_)
        throw std::runtime_error(str(format("Failed to load bam index for %1%") % path));

    if (bam_parse_region(header(), region, &region_.tid, &region_.beg, &region_.end) < 0) {
        throw std::runtime_error(str(format(
            "Failed to parse bam region '%1%' in file %2%. ")
            % region % path));
    }

    iter_ = bam_iter_query(index_, region_.tid, region_.beg, region_.end);
}

inline
Region const* RegionLimitedBamReader::region() const {
    return &region_;
}

inline
RegionLimitedBamReader::~RegionLimitedBamReader() {
    bam_iter_destroy(iter_);
    bam_index_destroy(index_);
    index_ = 0;
}

inline
bool RegionLimitedBamReader::takeImpl(bam1_t* entry) {
    return bam_iter_read(BamReader::fp_->x.bam, iter_, entry) > 0;
}
