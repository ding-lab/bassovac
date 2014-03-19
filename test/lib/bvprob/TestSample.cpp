#include "bvprob/Sample.hpp"
#include "bvprob/PBin.hpp"

#include <gtest/gtest.h>

using namespace std;

TEST(TestSample, phi) {
    double variantFrequency = 0.001;
    Sample s;
    s.setValues(10, 5, variantFrequency, 1.0, 1.0, 0, 0, 0);

    AlleleType alleles[2] = { REF, VAR };
    ASSERT_EQ(variantFrequency, s.piecewisePhi(alleles));
    alleles[1] = REF; // REF, REF
    ASSERT_EQ(0.0, s.piecewisePhi(alleles));
    alleles[0] = alleles[1] = VAR;
    ASSERT_EQ(1.0, s.piecewisePhi(alleles));
}
