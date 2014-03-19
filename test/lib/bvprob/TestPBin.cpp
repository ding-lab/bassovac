#include "bvprob/PBin.hpp"
#include "utility/Lut.hpp"

#include <gtest/gtest.h>

#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <ctime>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <vector>

using namespace std;

namespace {
    ostream& operator<<(ostream& s, const vector<double>& v) {
        s.flags(ios::fixed);
        s.fill(' ');
        s.width(9);
        s.precision(7);
        ostream_iterator<double> iter(s, ", ");
        copy(v.begin(), v.end(), iter);
        return s;
    }
}

TEST(TestPBin, binPValues2) {
    Lut::init();

    uint8_t probs[] = {
        2,2,2,2, // bin 1
        14,14,   // bin 2
        23,23,
        26,27,28,31,31,32,33,33,
        37,37,37,37,37,37,37,37,38,39,39,39
    };
    struct {
        uint32_t begin;
        uint32_t end;
    } expectedRanges[] = {
        { 0, 4 },
        { 4, 28 }
    };

    uint32_t expectedSizes[] = {
        4, 24
    };

    unsigned nProbs = sizeof(probs)/sizeof(probs[0]);

    size_t maxBins = 2;
    vector<PBin> bins(maxBins);
    uint32_t nBins = binPValues(probs, nProbs, maxBins, bins);
    ASSERT_EQ(maxBins, nBins);
    for (unsigned i = 0; i < maxBins; ++i) {
        vector<double> expected(
            probs+expectedRanges[i].begin,
            probs+expectedRanges[i].end
        );
        ASSERT_EQ(expectedSizes[i], bins[i].size)
            << "in bin " << i;
    }
}

TEST(TestPBin, binPValues) {
    Lut::init();

    uint8_t probs[] = {
        2,2,2,2, // bin 1
        14,14,   // bin 2
        23,23,   // bin 3 ...
        26,27,28,31,31,32,33,33,
        37,37,37,37,37,37,37,37,38,39,39,39
    };
    struct {
        uint32_t begin;
        uint32_t end;
    } expectedRanges[] = {
        { 0, 4 },
        { 4, 6 },
        { 6, 8 },
        { 8, 16 },
        { 16, 28 },
    };

    uint32_t expectedSizes[] = {
        4, 2, 2, 8, 12
    };

    unsigned nProbs = sizeof(probs)/sizeof(probs[0]);

    vector<PBin> bins(5);
    uint32_t nBins = binPValues(probs, nProbs, 5, bins);
    ASSERT_EQ(5u, nBins);
    for (unsigned i = 0; i < 5; ++i) {
        vector<double> expected(
            probs+expectedRanges[i].begin,
            probs+expectedRanges[i].end
        );
        ASSERT_EQ(expectedSizes[i], bins[i].size)
            << "in bin " << i;
    }
}
