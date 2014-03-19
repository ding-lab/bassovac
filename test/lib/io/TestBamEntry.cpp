#include "io/BamEntry.hpp"
#include "io/BamReader.hpp"
#include "utility/TempFile.hpp"

#include <bam.h>

#include <gtest/gtest.h>

#include <cstdint>
#include <cstring>
#include <iostream>
#include <string>
#include <vector>

namespace {
    std::string const SAM_DATA =
        "@SQ\tSN:1\tLN:247249719\n"
        "@SQ\tSN:2\tLN:247249719\n"
        "@SQ\tSN:3\tLN:247249719\n"
        "READ1\t73\t1\t101\t10\t10M\t=\t118985\t0\tACGTNACGTN\t0123456789\n"
        "READ2\t73\t2\t201\t10\t2M1I3M2D4M\t=\t118985\t0\tACGTNACGTN\t0123456789\n"
        "READ3\t73\t3\t201\t10\t2S3M1I4M\t=\t118985\t0\tACGTNACGTN\t0123456789\n"
        ;
}

std::ostream& operator<<(std::ostream& s, BamEntry::PileupData const& x) {
    s << "PileupData(tid=" << x.tid << ", base=" << x.base
        << ", quality=" << x.quality << ")";
    return s;
}


class TestBamEntry : public ::testing::Test {
public:
    void SetUp() {
        TempDir tmpdir;
        auto samFile = tmpdir.tempFile(SAM_DATA);
        BamReader reader(samFile->path());

        while (BamEntry* e = reader.take()) {
            entries.push_back(e);
        }
    }

    void TearDown() {
        for (auto i = entries.begin(); i != entries.end(); ++i) {
            delete *i;
        }
    }

    std::vector<BamEntry*> entries;
};

TEST_F(TestBamEntry, parseCigarAllMatch) {
    auto const& e = *entries[0]; // cigar = 10M

    BamEntry::PileupData expected;
    expected.tid = 0;
    expected.base = 1;
    expected.quality = int('0') - 33;

    int bases[] = {1, 2, 4, 8, 15};

    BamEntry::PileupData x;
    for (int i = 0; i < 10; ++i) {
        expected.base = bases[i % 5];
        EXPECT_TRUE(e.resolveCigar(100 + i, x));
        EXPECT_EQ(expected, x) << "at base index " << i;
        ++expected.readOffset;
        ++expected.quality;
    }
}

TEST_F(TestBamEntry, parseCigarWithInsertion) {
    auto const& e = *entries[1]; // cigar = 2M1I3M2D4M

    BamEntry::PileupData x;
    EXPECT_TRUE(e.resolveCigar(200, x));
    EXPECT_EQ(1, x.tid);
    EXPECT_EQ(1, x.base);

    EXPECT_TRUE(e.resolveCigar(201, x));
    EXPECT_EQ(1, x.tid);
    EXPECT_EQ(2, x.base);

    // Now we're skipping ahead over the insertion, and we're on the 3M
    EXPECT_TRUE(e.resolveCigar(202, x));
    EXPECT_EQ(1, x.tid);
    EXPECT_EQ(8, x.base);

    EXPECT_TRUE(e.resolveCigar(203, x));
    EXPECT_EQ(1, x.tid);
    EXPECT_EQ(15, x.base);

    EXPECT_TRUE(e.resolveCigar(204, x));
    EXPECT_EQ(1, x.tid);
    EXPECT_EQ(1, x.base);

    // 2D
    EXPECT_FALSE(e.resolveCigar(205, x));
    EXPECT_FALSE(e.resolveCigar(206, x));

    // 4M
    int expectedBases[] = {2, 4, 8, 15};
    for (int i = 0; i < 4; ++i) {
        EXPECT_TRUE(e.resolveCigar(207 + i, x));
        EXPECT_EQ(1, x.tid);
        EXPECT_EQ(expectedBases[i], x.base);
    }
}

TEST_F(TestBamEntry, resolveCigarSkippingCigarEntries) {
    auto const& e = *entries[1]; // cigar = 2M1I3M2D4M

    BamEntry::PileupData x;
    EXPECT_TRUE(e.resolveCigar(200, x));
    EXPECT_EQ(1, x.tid);
    EXPECT_EQ(1, x.base);

    // 4M
    int expectedBases[] = {2, 4, 8, 15};
    for (int i = 0; i < 4; ++i) {
        EXPECT_TRUE(e.resolveCigar(207 + i, x));
        EXPECT_EQ(1, x.tid);
        EXPECT_EQ(expectedBases[i], x.base) << "Iter " << i;
    }
}

TEST_F(TestBamEntry, resolveCigarSoftClipLeft) {
    auto const& e = *entries[2]; // cigar = 2S3M1I2D4M

    BamEntry::PileupData x;
    EXPECT_TRUE(e.resolveCigar(200, x));
    EXPECT_EQ(2, x.tid);
    EXPECT_EQ(4, x.base);

    EXPECT_TRUE(e.resolveCigar(201, x));
    EXPECT_EQ(8, x.base);

    EXPECT_TRUE(e.resolveCigar(202, x));
    EXPECT_EQ(15, x.base);

    EXPECT_TRUE(e.resolveCigar(203, x));
    EXPECT_EQ(2, x.base);

    EXPECT_TRUE(e.resolveCigar(204, x));
    EXPECT_EQ(4, x.base);
}
