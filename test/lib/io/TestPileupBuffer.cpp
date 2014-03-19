#include "io/PileupBuffer.hpp"
#include "io/Pileup.hpp"
#include "io/BamReader.hpp"
#include "utility/TempFile.hpp"

#include <gtest/gtest.h>
#include <stdexcept>
#include <string>
#include <memory>

using namespace std;

namespace {
    typedef shared_ptr<Pileup> PileupPtr;
}

class TestPileupBuffer : public testing::Test {
public:
    TempDir tmpdir;
};

TEST_F(TestPileupBuffer, chrom) {
    const string Sam =
        "@SQ\tSN:1\tLN:247249719\n"
        "@SQ\tSN:2\tLN:247249719\n"
        "@SQ\tSN:3\tLN:247249719\n"
        "READ1\t73\t2\t118985\t10\t36M\t=\t118985\t0\tATAAAAATCTATCATTTCTCCTTCCAGTTTTTTTTT\t2<<<9<<<<<7<&7<<7<7<<9<<<**82%.31,4,\n"
        "READ1\t73\t1\t118985\t10\t36M\t=\t118985\t0\tATAAAAATCTATCATTTCTCCTTCCAGTTTTTTTTT\t2<<<9<<<<<7<&7<<7<7<<9<<<**82%.31,4,\n"
        "READ1\t73\t3\t118985\t10\t36M\t=\t118985\t0\tATAAAAATCTATCATTTCTCCTTCCAGTTTTTTTTT\t2<<<9<<<<<7<&7<<7<7<<9<<<**82%.31,4,\n"
        "READ1\t73\t2\t118985\t10\t36M\t=\t118985\t0\tATAAAAATCTATCATTTCTCCTTCCAGTTTTTTTTT\t2<<<9<<<<<7<&7<<7<7<<9<<<**82%.31,4,\n"
        ;

    auto samFile(tmpdir.tempFile(Sam));
    BamReader reader(samFile->path());

    PileupBuffer pb;
    ASSERT_TRUE(pb.push(reader.take()));

    std::unique_ptr<BamEntry> entry(reader.take());
    ASSERT_FALSE(pb.push(entry.get()));
    entry.reset(reader.take());
    ASSERT_FALSE(pb.push(entry.get()));

    ASSERT_TRUE(pb.push(reader.take()));
}

TEST_F(TestPileupBuffer, disjoint) {
    const string Sam =
        "@SQ\tSN:1\tLN:247249719\n"
        "READ1\t133\t1\t1\t0\t4M\t=\t118985\t0\tAATT\t20<<\n"
        "READ2\t73\t1\t2\t10\t1M\t=\t118985\t0\tA\t2\n"
        "READ3\t73\t1\t5\t10\t4M\t=\t118985\t0\tATAA\t2<<<\n"
        ;

    auto samFile = tmpdir.tempFile(Sam);
    BamReader reader(samFile->path());

    PileupBuffer pb;
    ASSERT_TRUE(pb.push(reader.take()));
    ASSERT_TRUE(pb.push(reader.take()));
    std::unique_ptr<BamEntry> entry(reader.take());
    ASSERT_FALSE(pb.push(entry.get()));
}

TEST_F(TestPileupBuffer, clear) {

    const string Sam =
        "@SQ\tSN:1\tLN:247249719\n"
        "@SQ\tSN:2\tLN:247249719\n"
        "READ1\t73\t2\t118985\t10\t36M\t=\t118985\t0\tATAAAAATCTATCATTTCTCCTTCCAGTTTTTTTTT\t<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\n"
        "READ1\t73\t1\t118985\t10\t36M\t=\t118985\t0\tATAAAAATCTATCATTTCTCCTTCCAGTTTTTTTTT\t>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n"
        ;

    auto samFile = tmpdir.tempFile(Sam);
    BamReader reader(samFile->path());

    PileupBuffer pb;
    ASSERT_TRUE(pb.push(reader.take()));
    std::unique_ptr<BamEntry> entry(reader.take());
    ASSERT_FALSE(pb.push(entry.get()));

    pb.clear();
    ASSERT_TRUE(pb.empty());
    ASSERT_EQ(0, pb.start());
    ASSERT_EQ(0u, pb.end());
}

TEST_F(TestPileupBuffer, clearBefore) {
    const string Sam =
        "@SQ\tSN:1\tLN:247249719\n"
        "READ1\t133\t1\t1\t0\t4M\t=\t118985\t0\tAATT\t<<<<\n"
        "READ2\t73\t1\t3\t10\t8M\t=\t118985\t0\tAAAAAAAA\t<<<<<<<<\n"
        "READ3\t73\t1\t4\t10\t8M\t=\t118985\t0\tCCCCCCCC\t>>>>>>>>\n"
        ;

    auto samFile = tmpdir.tempFile(Sam);
    BamReader reader(samFile->path());

    PileupBuffer pb;
    ASSERT_TRUE(pb.push(reader.take()));
    ASSERT_TRUE(pb.push(reader.take()));
    ASSERT_TRUE(pb.push(reader.take()));

    pb.clearBefore(0, 5);
    // 0-based
    ASSERT_EQ(2u, pb.size());
    ASSERT_EQ(10u, pb.end()); // READ2 is now first, and its end is 10

    pb.clearBefore(1, 0);
    ASSERT_EQ(0u, pb.size());
}

TEST_F(TestPileupBuffer, cmpOverlap) {
    int32_t start = 0;
    uint32_t end = 0;

    const string Sam =
        "@SQ\tSN:1\tLN:247249719\n"
        "READ1\t133\t1\t10\t0\t4M\t=\t118985\t0\tAATT\t<<<<\n"
        "READ2\t73\t1\t13\t10\t8M\t=\t118985\t0\tAAAAAAAA\t<<<<<<<<\n"
        ;

    auto samFile = tmpdir.tempFile(Sam);
    BamReader reader(samFile->path());

    PileupBuffer pb1;
    PileupBuffer pb2;
    ASSERT_TRUE(pb1.push(reader.take()));
    ASSERT_TRUE(pb2.push(reader.take()));
    ASSERT_EQ(OVERLAP, pb1.cmp(pb2, &start, &end));
    // 0-based
    ASSERT_EQ(12, start);
    ASSERT_EQ(13u, end);
}

TEST_F(TestPileupBuffer, cmpAfter) {
    int32_t start;
    uint32_t end;

    const string Sam =
        "@SQ\tSN:1\tLN:247249719\n"
        "READ1\t133\t1\t9\t0\t4M\t=\t118985\t0\tAATT\t<<<<\n"
        "READ2\t73\t1\t1\t10\t8M\t=\t118985\t0\tAAAAAAAA\t<<<<<<<<\n"
        ;

    auto samFile = tmpdir.tempFile(Sam);
    BamReader reader(samFile->path());

    PileupBuffer pb1;
    PileupBuffer pb2;
    ASSERT_TRUE(pb1.push(reader.take()));
    ASSERT_TRUE(pb2.push(reader.take()));

    ASSERT_EQ(AFTER, pb1.cmp(pb2, &start, &end));
}

TEST_F(TestPileupBuffer, cmpBefore) {
    int32_t start;
    uint32_t end;

    const string Sam =
        "@SQ\tSN:1\tLN:247249719\n"
        "READ1\t133\t1\t1\t0\t4M\t=\t118985\t0\tAATT\t<<<<\n"
        "READ2\t73\t1\t5\t10\t8M\t=\t118985\t0\tAAAAAAAA\t<<<<<<<<\n"
        ;

    auto samFile = tmpdir.tempFile(Sam);
    BamReader reader(samFile->path());

    PileupBuffer pb1;
    PileupBuffer pb2;
    ASSERT_TRUE(pb1.push(reader.take()));
    ASSERT_TRUE(pb2.push(reader.take()));

    ASSERT_EQ(BEFORE, pb1.cmp(pb2, &start, &end));
}

TEST_F(TestPileupBuffer, testPileup) {
    // A note on the quality values (last column):
    // numeric (phred) quality will equal ord(c)-33 where ord(c) is the ascii value of a
    // particular character in the string (e.g., B=ord('A')-33=65-33=32.
    const string Sam =
        "@SQ\tSN:1\tLN:247249719\n"
        "READ1\t83\t1\t2\t0\t3M\t=\t118985\t0\tCAT\tABC\n"
        "READ2\t83\t1\t3\t10\t3M\t=\t118985\t0\tCAT\tABC\n"
        "READ3\t83\t1\t4\t10\t3M\t=\t118985\t0\tCAT\tABC\n"
        ;

    auto samFile = tmpdir.tempFile(Sam);
    BamReader reader(samFile->path());

    PileupBuffer pb;
    ASSERT_TRUE(pb.push(reader.take()));
    ASSERT_TRUE(pb.push(reader.take()));
    ASSERT_TRUE(pb.push(reader.take()));

    PileupPtr reads(pb.pileup(0));
    ASSERT_TRUE(reads->empty());

    reads.reset(pb.pileup(1));
    ASSERT_EQ(1u, reads->size());
//    ASSERT_EQ(string("READ1"), (*reads)[0]->name());
    ASSERT_EQ('C', bam_nt16_rev_table[(*reads)[0].base]);
    ASSERT_EQ(32, (*reads)[0].quality);

    reads.reset(pb.pileup(2));
    ASSERT_EQ(2u, reads->size());
/*
    ASSERT_EQ(string("READ1"), (*reads)[0]->name());
    ASSERT_EQ(string("READ2"), (*reads)[1]->name());
*/
    ASSERT_EQ('A', bam_nt16_rev_table[(*reads)[0].base]);
    ASSERT_EQ('C', bam_nt16_rev_table[(*reads)[1].base]);
    ASSERT_EQ(33, (*reads)[0].quality);
    ASSERT_EQ(32, (*reads)[1].quality);

    reads.reset(pb.pileup(3));
    ASSERT_EQ(3u, reads->size());
/*
    ASSERT_EQ(string("READ1"), (*reads)[0]->name());
    ASSERT_EQ(string("READ2"), (*reads)[1]->name());
    ASSERT_EQ(string("READ3"), (*reads)[2]->name());
*/
    ASSERT_EQ('T', bam_nt16_rev_table[(*reads)[0].base]);
    ASSERT_EQ('A', bam_nt16_rev_table[(*reads)[1].base]);
    ASSERT_EQ('C', bam_nt16_rev_table[(*reads)[2].base]);
    ASSERT_EQ(34, (*reads)[0].quality);
    ASSERT_EQ(33, (*reads)[1].quality);
    ASSERT_EQ(32, (*reads)[2].quality);
    int baseCounts[4] = {0};
    reads->baseCounts(baseCounts);
    ASSERT_EQ(1, baseCounts[0]);
    ASSERT_EQ(1, baseCounts[1]);
    ASSERT_EQ(0, baseCounts[2]);
    ASSERT_EQ(1, baseCounts[3]);

    reads.reset(pb.pileup(4));
    ASSERT_EQ(2u, reads->size());
/*
    ASSERT_EQ(string("READ2"), (*reads)[0]->name());
    ASSERT_EQ(string("READ3"), (*reads)[1]->name());
*/
    ASSERT_EQ('T', bam_nt16_rev_table[(*reads)[0].base]);
    ASSERT_EQ('A', bam_nt16_rev_table[(*reads)[1].base]);
    ASSERT_EQ(34, (*reads)[0].quality);
    ASSERT_EQ(33, (*reads)[1].quality);

    reads.reset(pb.pileup(5));
    ASSERT_EQ(1u, reads->size());
//    ASSERT_EQ(string("READ3"), (*reads)[0]->name());
    ASSERT_EQ('T', bam_nt16_rev_table[(*reads)[0].base]);
    ASSERT_EQ(34, (*reads)[0].quality);

    reads.reset(pb.pileup(6));
    ASSERT_TRUE(reads->empty());
}

TEST_F(TestPileupBuffer, testDeletion) {
    const string Sam =
        "@SQ\tSN:1\tLN:247249719\n"
        "READ1\t83\t1\t2\t0\t3M\t=\t118985\t0\tCAT\tABC\n"
        "READ2\t83\t1\t2\t10\t1M1D2M\t=\t118985\t0\tCAT\tABC\n"
        ;

    auto samFile = tmpdir.tempFile(Sam);
    BamReader reader(samFile->path());

    PileupBuffer pb;
    ASSERT_TRUE(pb.push(reader.take()));
    ASSERT_TRUE(pb.push(reader.take()));

    PileupPtr reads(pb.pileup(0));
    ASSERT_TRUE(reads->empty());

    reads.reset(pb.pileup(1));
    ASSERT_EQ(2u, reads->size());

    reads.reset(pb.pileup(2));
    ASSERT_EQ(1u, reads->size());

    reads.reset(pb.pileup(3));
    ASSERT_EQ(2u, reads->size());

    reads.reset(pb.pileup(4));
    ASSERT_EQ(1u, reads->size());
}
