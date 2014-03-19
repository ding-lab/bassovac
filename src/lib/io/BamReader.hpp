#pragma once

#include "BamReaderBase.hpp"

#include <boost/format.hpp>
#include <functional>
#include <stdexcept>
#include <string>

class BamReader : public BamReaderBase {
public:
    explicit BamReader(std::string const& path);
    ~BamReader();

    bam_header_t* header() const;
    std::string const& path() const;

protected:
    virtual bool takeImpl(bam1_t* entry);

protected:
    std::string path_;
    samfile_t* fp_;
};

inline
BamReader::BamReader(std::string const& path)
    : path_(path)
    , fp_(0)
{
    using boost::format;

    std::string mode = "r";
    if (path_.find(".bam") != std::string::npos)
        mode += 'b';

    fp_ = samopen(path_.c_str(), mode.c_str(), 0);
    if (!fp_ || !fp_->header) {
        throw std::runtime_error(str(format("Failed to open samfile %1%") % path));
    }

    if (!fp_->x.bam) {
        throw std::runtime_error(str(format("%1% is not a valid bam file") % path));
    }
}

inline
BamReader::~BamReader() {
    samclose(fp_);
    fp_ = 0;
}

inline
bool BamReader::takeImpl(bam1_t* entry) {
    return samread(fp_, entry) > 0;
}

inline
bam_header_t* BamReader::header() const {
    return fp_->header;
}

inline
std::string const& BamReader::path() const {
    return path_;
}
