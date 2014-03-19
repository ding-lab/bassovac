#include "utility/Binomial.hpp"
#include "utility/Lut.hpp"

#include <gtest/gtest.h>

#include <cstdint>
#include <cstdlib>
#include <ctime>

namespace {
    // These values are generated by the R script pdfConvolve2 that lives
    // next to this file.
    struct Params {
        double p1, p2;
        size_t n1, n2, k;
        double result;
    };

    Params params[] = {
        //p1            p2            n1   n2   k  result
        { 0.4019479589, 0.839249833,   5,   8,  6, 5.310741e-02 },
        { 0.5782658223, 0.356853656,  20,  20, 20, 1.177496e-01 },
        { 0.2780696503, 0.696222218,  16,  15, 15, 1.570055e-01 },
        { 0.0005520571, 0.541839660,  13,  11, 12, 8.688104e-06 },
        { 0.9164467736, 0.412507781,   7,   2,  4, 5.618983e-03 },
        { 0.2702691413, 0.841355638,  16,   4, 10, 9.696953e-02 },
        { 0.1065676345, 0.992481831,  11,   7,  9, 2.189948e-01 },
        { 0.0768027040, 0.205165485,  13,  10, 11, 3.639990e-05 },
        { 0.6778377453, 0.581781673,  15,  11, 13, 5.574166e-02 },
        { 0.7582510132, 0.001995862,  14,  13, 13, 9.526076e-02 },
    };

    size_t const n = sizeof(params) / sizeof(params[0]);
}

class TestBinomial : public ::testing::Test {
public:
    void SetUp() {
        Lut::init();
    }

};

TEST_F(TestBinomial, pdfConvolve2WorkedExample) {
    // Suppose X ~ Binom(1/2, 3), Y ~ Binom(1/4, 5), Z = X + Y
    //  x | P(X=x)
    // ---|-------
    //  0 | 1 * (1/2)^0 (1/2)^3 = 1/8
    //  1 | 3 * (1/2)^2 (1/2)^1 = 3/8
    //  2 | 3 * (1/2)^1 (1/2)^2 = 3/8
    //  3 | 1 * (1/2)^0 (1/2)^3 = 1/8
    // -----------
    //  y | P(Y=y)
    // ---|-------
    //  0 |  1 * (1/4)^0 (3/4)^5 = 243/1024
    //  1 |  5 * (1/4)^1 (3/4)^4 = 405/1024
    //  2 | 10 * (1/4)^2 (3/4)^3 = 270/1024
    //  3 | 10 * (1/4)^3 (3/4)^2 =  90/1024
    //  4 |  5 * (1/4)^4 (3/4)^1 =  15/1024
    //  5 |  1 * (1/4)^5 (3/4)^0 =   1/1024
    //
    // Then P(Z=3)
    //  = P(X=0,Y=3) + P(X=1,Y=2) + P(X=2,Y=1) + P(X=3,Y=0)
    //  = (1/8)(90/1024) + (3/8)(270/1024) + (3/8)(405/1024) + (1/8)(243/1024)
    //  = 45/4096 + 405/4096 + 1215/8192 + 243/8192
    //  = 2358/8192

    // Let's see what the function thinks...
    double p1 = 0.5;
    double p2 = 0.25;
    uint32_t n1 = 3;
    uint32_t n2 = 5;

    double p = Binomial::pdfConvolve2(p1, p2, n1, n2, 3);
    EXPECT_DOUBLE_EQ(2358.0 / 8192.0, p);
}

TEST_F(TestBinomial, pdfConvolve2CompareToR) {

    for (size_t i = 0; i < n; ++i) {
        auto const& p = params[i];
        double result = Binomial::pdfConvolve2(p.p1, p.p2, p.n1, p.n2, p.k);
        EXPECT_NEAR(p.result, result, 1e-7)
            << "p1=" << p.p1 << "; "
            << "p2=" << p.p2 << "; "
            << "n1=" << p.n1 << "; "
            << "n2=" << p.n2 << "; "
            << "k=" << p.k
            ;
    }
}