#include "SamConvert.hpp"
#include "BamEntry.hpp"
#include "BamReader.hpp"

#include <bam.h>

#include <boost/format.hpp>

#include <stdexcept>

using boost::format;

void samToIndexedBam(std::string const& samPath, std::string const& bamPath) {
    BamReader in(samPath);
    auto fp = samopen(bamPath.c_str(), "wb", in.header());
    if (!fp) {
        throw std::runtime_error(str(format(
            "Failed to open bam file %1% for writing!"
            ) % bamPath));
    }
    while (auto e = in.take()) {
        int rv = samwrite(fp, e->rawData());
        delete e;
        if (rv < 0) {
            throw std::runtime_error(str(format(
                "Failed to write bam entry to %1%"
                ) % bamPath));
        }
    }
    samclose(fp);
    bam_index_build(bamPath.c_str());
}
