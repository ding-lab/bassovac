#include "io/BamIntersector.hpp"
#include "io/BamReader.hpp"
#include "io/RegionLimitedBamReader.hpp"
#include "io/Pileup.hpp"
#include "io/SamConvert.hpp"
#include "utility/TempFile.hpp"

#include <gtest/gtest.h>

#include <deque>
#include <functional>
#include <map>
#include <set>

using namespace std;
using namespace std::placeholders;

namespace {
    const string normalSam =
        "@SQ\tSN:1\tLN:247249719\n"
        "READ1\t83\t1\t1\t60\t10M\t=\t118985\t0\tAAAAAAAAAA\t<<<<<<<<<<\n"
        "READ2\t83\t1\t9\t60\t4M\t=\t118985\t0\tGGGG\t<<<<\n"
        "READ3\t83\t1\t11\t60\t10M\t=\t118985\t0\tCCCCCCCCCC\t>>>>>>>>>>\n"
        "READ4\t83\t1\t19\t60\t6M\t=\t118985\t0\tTTTTTT\t>>>>>>\n"
        "READ5\t83\t1\t22\t60\t10M\t=\t118985\t0\tAAAAAAAAAA\t>>>>>>>>>>\n"
        ;

    const string tumorSam =
        "@SQ\tSN:1\tLN:247249719\n"
        "READA\t83\t1\t5\t60\t6M\t=\t118985\t0\tMMMMMM\t<<<<<<\n"
        "READB\t83\t1\t6\t60\t6M\t=\t118985\t0\tSSSSSS\t<<<<<<\n"
        "READC\t83\t1\t12\t60\t4M\t=\t118985\t0\tWWWW\t>>>>\n"
        "READD\t83\t1\t12\t60\t5M\t=\t118985\t0\tDDDDD\t<<<<<\n"
        "READE\t83\t1\t19\t60\t10M\t=\t118985\t0\tNNNNNNNNNN\t>>>>>>>>>>\n"
        "READF\t83\t1\t21\t60\t10M\t=\t118985\t0\tKKKKKKKKKK\t<<<<<<<<<<\n"
        ;


    struct ReadCounts {
        ReadCounts() : normalCount(0), tumorCount(0) {}
        ReadCounts(uint32_t n, uint32_t t) : normalCount(n), tumorCount(t) {}
        uint32_t normalCount;
        uint32_t tumorCount;
    };

    struct Collector {
        typedef pair<uint32_t, Pileup> PosReads;

        void collect(int32_t pos, const Pileup& normal, const Pileup& tumor) {
            results[pos] = ReadCounts(normal.size(), tumor.size());
        }

        map<int32_t, ReadCounts> results;
    };
}

class TestBamIntersector : public ::testing::Test {
public:
    void SetUp() {
        auto normalFile(tmpdir.tempFile(normalSam));
        auto tumorFile(tmpdir.tempFile(tumorSam));

        normalBamPath = normalFile->path() + ".bam";
        tumorBamPath = tumorFile->path() + ".bam";
        samToIndexedBam(normalFile->path(), normalBamPath);
        samToIndexedBam(tumorFile->path(), tumorBamPath);
    }

protected:
    TempDir tmpdir;
    std::string normalBamPath;
    std::string tumorBamPath;
};

TEST_F(TestBamIntersector, intersect) {
    BamReader normalReader(normalBamPath);
    BamReader tumorReader(tumorBamPath);
    Collector collector;

    BamIntersector intersector(normalReader, tumorReader,
        std::bind(&Collector::collect, &collector, _1, _2, _3));
    intersector.run();

// Test reads look like:
//   0123456789012345678901234567890
// A [--------][--------] [--------]
// A         [--]      [----]
//
// B     [----] [--]   [--------]
// B      [----][---]    [--------]

    // these 'expected' arrays count the number of expected reads at each position
    uint32_t expectedNormal[] = {
        0, 0, 0, 0, 1,
        1, 1, 1, 2, 2,
        2, 2, 1, 1, 1,
        1, 1, 1, 2, 2,
        1, 2, 2, 2, 1,
        1, 1, 1, 1, 1
    };

    uint32_t expectedTumor[] = {
        0, 0, 0, 0, 1,
        2, 2, 2, 2, 2,
        1, 2, 2, 2, 2,
        1, 0, 0, 1, 1,
        2, 2, 2, 2, 2,
        2, 2, 2, 1, 1
    };

    for (auto iter = collector.results.begin(); iter != collector.results.end(); ++iter) {
        ASSERT_EQ(expectedNormal[iter->first], iter->second.normalCount) << "position " << iter->first;
        ASSERT_EQ(expectedTumor[iter->first], iter->second.tumorCount) << "position " << iter->first;
    }
}

TEST_F(TestBamIntersector, intersectRegionLimited) {
    std::string region("1:9-11");
    // FIXME: go make RLBamReader take std::string for region
    RegionLimitedBamReader normalReader(normalBamPath, region.c_str());
    RegionLimitedBamReader tumorReader(tumorBamPath, region.c_str());

    Collector collector;

    BamIntersector intersector(normalReader, tumorReader,
        std::bind(&Collector::collect, &collector, _1, _2, _3));
    intersector.run();

    // Test reads look like, well, the same as they do in the previous test.
    // Now we just want to see that the positions we observe are actually
    // limited by the region. The bam reader is going to return all reads that
    // intersect the region, so the intersector has to do some limiting of its
    // own.
    std::set<uint32_t> observedPositions;
    for (auto iter = collector.results.begin(); iter != collector.results.end(); ++iter) {
        observedPositions.insert(iter->first);
    }
    std::set<uint32_t> expected;
    expected.insert(8);
    expected.insert(9);
    expected.insert(10);

    EXPECT_EQ(expected, observedPositions);
}


