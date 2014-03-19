#include "io/BamEntry.hpp"
#include "io/BamFilter.hpp"
#include "io/BamReader.hpp"
#include "io/RegionLimitedBamReader.hpp"
#include "io/SamConvert.hpp"
#include "utility/TempFile.hpp"

#include <boost/filesystem.hpp>
#include <gtest/gtest.h>
#include <cstdlib>
#include <fstream>
#include <stdexcept>
#include <string>

namespace bfs = boost::filesystem;
using namespace std;

namespace {
    const string samSource =
        "@SQ\tSN:1\tLN:247249719\n"
        "@SQ\tSN:2\tLN:247249719\n"
        "READ1\t73\t1\t118985\t10\t36M\t=\t118985\t0\tATAAAAATCTATCATTTCTCCTTCCAGTTTTTTTTT\t2<<<9<<<<<7<&7<<7<7<<9<<<**82%.31,4,\n"
        "READ2\t133\t1\t218985\t0\t*\t=\t118985\t0\tAATTCCATCTCTGCGTCTTTCCTCCTTCTTCTCTTC\t20<<<<2<2<<<4&<3*0<3<<<<<6*<*&*0+0%'\n"
        "READ3\t73\t2\t118985\t20\t36M\t=\t118985\t0\tATAAAAATCTATCATTTCTCCTTCCAGTTTTTTTTT\t2<<<9<<<<<7<&7<<7<7<<9<<<**82%.31,4,\n"
        ;
}


class TestBamReader : public ::testing::Test {
public:
    void SetUp() {
        samFile = tmpdir.tempFile(samSource);
    }

    std::string makeBamFromSam(std::string const& path) {
        std::string output = path + ".bam";
        samToIndexedBam(path, output);
        return output;
    }

protected:
    TempDir tmpdir;
    std::unique_ptr<TempFile> samFile;
    std::unique_ptr<TempFile> bamFile;
};

TEST_F(TestBamReader, read) {
    BamReader reader(samFile->path());

    // no rhyme or region to it.
    ASSERT_FALSE(reader.region());

    BamEntry* e = reader.peek();
    ASSERT_TRUE(e);
    ASSERT_STREQ("READ1", e->name());
    ASSERT_EQ(118984, e->start());
    ASSERT_EQ(119020, e->end());

    ASSERT_STREQ("1", reader.targetName(e->tid()));

    e = reader.peek();
    ASSERT_STREQ("READ1", e->name());

    e = reader.take();
    ASSERT_STREQ("READ1", e->name());
    delete e;

    e = reader.peek();
    ASSERT_STREQ("READ2", e->name());

    e = reader.take();
    ASSERT_STREQ("READ2", e->name());
    delete e;

    e = reader.take();
    ASSERT_STREQ("READ3", e->name());
    delete e;

    e = reader.peek();
    ASSERT_FALSE(e);
    ASSERT_FALSE(reader.take());
    ASSERT_FALSE(reader.take());
}

TEST_F(TestBamReader, filtered10) {
    BamFilter filter(BAM_DEF_MASK, 10);
    BamReader reader(samFile->path());
    reader.setFilter(&filter);

    BamEntry* e = reader.peek();
    ASSERT_TRUE(e);
    ASSERT_STREQ("READ1", e->name());

    e = reader.peek();
    ASSERT_STREQ("READ1", e->name());

    e = reader.take();
    ASSERT_STREQ("READ1", e->name());
    delete e;

    // read 2 is filtered
    e = reader.peek();
    ASSERT_STREQ("READ3", e->name());

    e = reader.take();
    ASSERT_STREQ("READ3", e->name());
    delete e;

    e = reader.peek();
    ASSERT_FALSE(e);
    ASSERT_FALSE(reader.take());
    ASSERT_FALSE(reader.take());
}

TEST_F(TestBamReader, filtered11) {
    BamFilter filter(BAM_DEF_MASK, 11);
    BamReader reader(samFile->path());
    reader.setFilter(&filter);

    // reads 1 & 2 are filtered
    auto e = reader.peek();
    ASSERT_STREQ("READ3", e->name());

    e = reader.take();
    ASSERT_STREQ("READ3", e->name());
    delete e;

    e = reader.peek();
    ASSERT_FALSE(e);
    ASSERT_FALSE(reader.take());
    ASSERT_FALSE(reader.take());
}

TEST_F(TestBamReader, regionLimitedAndFiltered) {
    auto bamFile = makeBamFromSam(samFile->path());

    // limit to all of chr 1
    std::unique_ptr<RegionLimitedBamReader> reader(new RegionLimitedBamReader(bamFile, "1"));
    EXPECT_EQ(0, reader->region()->tid);
    EXPECT_EQ(0, reader->region()->beg);

    // I think the actual value of 'end' is based on the max size of a sequence
    // that bam can index (2^29). What we'll do instead is make sure it's
    // >= the sequence length we specified in the header (247249719)
    EXPECT_GE(reader->region()->end, 247249719);

    auto e = reader->take();
    EXPECT_STREQ("READ1", e->name());
    delete e;
    e = reader->take();
    EXPECT_STREQ("READ2", e->name());
    delete e;
    ASSERT_FALSE(reader->take());
    ASSERT_FALSE(reader->peek());

    // limit to all of chr 1 with filter excluding READ2
    reader.reset(new RegionLimitedBamReader(bamFile, "1"));
    BamFilter filter(BAM_DEF_MASK, 1);
    reader->setFilter(&filter);
    e = reader->take();
    EXPECT_STREQ("READ1", e->name());
    delete e;
    ASSERT_FALSE(reader->take());
    ASSERT_FALSE(reader->peek());


    // limit to all of chr 2
    reader.reset(new RegionLimitedBamReader(bamFile, "2"));
    EXPECT_EQ(1, reader->region()->tid);
    EXPECT_EQ(0, reader->region()->beg);
    EXPECT_GE(reader->region()->end, 247249719);

    e = reader->take();
    EXPECT_STREQ("READ3", e->name());
    delete e;
    ASSERT_FALSE(reader->peek());
    ASSERT_FALSE(reader->take());

    // limit to part of chr 1
    // 118985
    reader.reset(new RegionLimitedBamReader(bamFile, "1:118985-118986"));
    EXPECT_EQ(0, reader->region()->tid);
    EXPECT_EQ(118984, reader->region()->beg); // zero based
    EXPECT_EQ(reader->region()->end, 118986);

    e = reader->take();
    EXPECT_STREQ("READ1", e->name());
    delete e;
    EXPECT_FALSE(reader->take());
}
