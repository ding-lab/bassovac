#include "io/CigarParser.hpp"

#include <gtest/gtest.h>

#include <bam.h>

#include <iostream>
#include <vector>

namespace {
    struct CigarBuilder {
        void add(int len, int op) {
            cigar.push_back(bam_cigar_gen(len, op));
        }

        std::vector<uint32_t> cigar;
    };
}

class TestCigarParser : public ::testing::Test {
public:
    CigarBuilder builder;
};

TEST_F(TestCigarParser, parseWithDeletion) {
    builder.add(2, BAM_CMATCH);
    builder.add(3, BAM_CDEL);
    builder.add(4, BAM_CMATCH);
    uint32_t const* cigar = builder.cigar.data();
    int n = builder.cigar.size();

    CigarParser p(cigar, n);

    for (int i = 0; i < 2; ++i) {
        EXPECT_EQ(0, p.currentOpIdx());
        EXPECT_EQ(BAM_CMATCH, p.currentOp());
        EXPECT_EQ(2, p.currentOpLen());
        EXPECT_EQ(0u, p.readPos());
        EXPECT_EQ(0u, p.refPos());

        EXPECT_EQ(0, p.getSnvPileupOffset(0));
        EXPECT_EQ(1, p.getSnvPileupOffset(1));
    }

    EXPECT_EQ(-1, p.getSnvPileupOffset(2));
    EXPECT_EQ(-1, p.getSnvPileupOffset(3));
    EXPECT_EQ(-1, p.getSnvPileupOffset(4));
    EXPECT_EQ( 2, p.getSnvPileupOffset(5));
    EXPECT_EQ( 3, p.getSnvPileupOffset(6));
    EXPECT_EQ( 4, p.getSnvPileupOffset(7));
    EXPECT_EQ( 5, p.getSnvPileupOffset(8));
}


TEST_F(TestCigarParser, parseWithInsertion) {
    builder.add(2, BAM_CMATCH);
    builder.add(3, BAM_CINS);
    builder.add(4, BAM_CMATCH);
    uint32_t const* cigar = builder.cigar.data();
    int n = builder.cigar.size();

    CigarParser p(cigar, n);

    for (int i = 0; i < 2; ++i) {
        EXPECT_EQ(0, p.currentOpIdx());
        EXPECT_EQ(BAM_CMATCH, p.currentOp());
        EXPECT_EQ(2, p.currentOpLen());
        EXPECT_EQ(0u, p.readPos());
        EXPECT_EQ(0u, p.refPos());

        EXPECT_EQ(0, p.getSnvPileupOffset(0));
        EXPECT_EQ(1, p.getSnvPileupOffset(1));
    }

    // read positions 2, 3, 4 are the insertion
    EXPECT_EQ( 5, p.getSnvPileupOffset(2));
    EXPECT_EQ( 6, p.getSnvPileupOffset(3));
    EXPECT_EQ( 7, p.getSnvPileupOffset(4));
    EXPECT_EQ( 8, p.getSnvPileupOffset(5));
}


TEST_F(TestCigarParser, parseWithSoftClipLeft) {
    builder.add(5, BAM_CSOFT_CLIP);
    builder.add(2, BAM_CMATCH);
    builder.add(3, BAM_CINS);
    builder.add(4, BAM_CMATCH);
    uint32_t const* cigar = builder.cigar.data();
    int n = builder.cigar.size();

    CigarParser p(cigar, n);

    EXPECT_EQ(0, p.currentOpIdx());
    EXPECT_EQ(BAM_CSOFT_CLIP, p.currentOp());
    EXPECT_EQ(5, p.currentOpLen());

    EXPECT_EQ(5, p.getSnvPileupOffset(0));
    EXPECT_EQ(6, p.getSnvPileupOffset(1));
    EXPECT_EQ(10, p.getSnvPileupOffset(2));
    EXPECT_EQ(11, p.getSnvPileupOffset(3));
    EXPECT_EQ(12, p.getSnvPileupOffset(4));
    EXPECT_EQ(13, p.getSnvPileupOffset(5));
}
