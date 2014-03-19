#include "utility/Lut.hpp"

#include <gtest/gtest.h>

#include <cmath>
#include <complex>

using namespace std;

TEST(TestLut, phred2p) {
    Lut::init();
    for (int i = 0; i < 256; ++i) {
        double p = Lut::phred2p(i);
        double pr = Lut::phred2p_reciprocal(i);
        ASSERT_DOUBLE_EQ(pow(10, i/-10.0), p);
        ASSERT_DOUBLE_EQ(pow(10, i/10.0), pr);
    }
}

TEST(TestLut, lgamma) {
    Lut::init(256);
    for (int i = 1; i < 256; ++i) {
        double p = Lut::lgamma(i);
        ASSERT_DOUBLE_EQ(std::lgamma(i), p);
    }
}

TEST(TestLut, rootsOfUnity) {
    Lut::init(100);
    EXPECT_EQ(std::complex<double>(1.0, 0.0), Lut::rootsOfUnity(1));
    for (unsigned i = 2; i <= 100; ++i) {
        auto expected = std::exp(std::complex<double>(0, 2.0 * M_PI / i));
        EXPECT_EQ(expected, Lut::rootsOfUnity(i)) << "at i=" << i;
    }
}
