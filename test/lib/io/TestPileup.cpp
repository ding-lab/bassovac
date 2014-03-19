#include "io/Pileup.hpp"

#include <gtest/gtest.h>

using namespace std;

// The baseCount arrays represent the counts of A,C,G,T
// The first param to variantAllele is the reference base (ACGT=powers of 2)
TEST(TestPileup, testVariantAllele) {
    {
        int baseCounts[4] = { 13, 0, 2, 3 };
        ASSERT_EQ(8, Pileup::variantAllele(1, baseCounts));
        ASSERT_EQ(1, Pileup::variantAllele(2, baseCounts));
        ASSERT_EQ(1, Pileup::variantAllele(4, baseCounts));
        ASSERT_EQ(1, Pileup::variantAllele(8, baseCounts));
    }

    {
        int baseCounts[4] = { 1, 4, 2, 3 };
        ASSERT_EQ(2, Pileup::variantAllele(1, baseCounts));
        ASSERT_EQ(8, Pileup::variantAllele(2, baseCounts));
        ASSERT_EQ(2, Pileup::variantAllele(4, baseCounts));
        ASSERT_EQ(2, Pileup::variantAllele(8, baseCounts));
    }

    {
        int baseCounts[4] = { 1, 1, 1, 1 };
        ASSERT_EQ(2, Pileup::variantAllele(1, baseCounts));
        ASSERT_EQ(1, Pileup::variantAllele(2, baseCounts));
        ASSERT_EQ(1, Pileup::variantAllele(4, baseCounts));
        ASSERT_EQ(1, Pileup::variantAllele(8, baseCounts));
    }

    {
        int baseCounts[4] = { 1, 0, 0, 0 };
        ASSERT_EQ(0, Pileup::variantAllele(1, baseCounts));
        ASSERT_EQ(1, Pileup::variantAllele(2, baseCounts));
        ASSERT_EQ(1, Pileup::variantAllele(4, baseCounts));
        ASSERT_EQ(1, Pileup::variantAllele(8, baseCounts));
    }

    {
        int baseCounts[4] = { 1, 1, 0, 0 };
        ASSERT_EQ(2, Pileup::variantAllele(1, baseCounts));
        ASSERT_EQ(1, Pileup::variantAllele(2, baseCounts));
        ASSERT_EQ(1, Pileup::variantAllele(4, baseCounts));
        ASSERT_EQ(1, Pileup::variantAllele(8, baseCounts));
        ASSERT_EQ(0, Pileup::variantAllele(3, baseCounts));
        ASSERT_EQ(1, Pileup::variantAllele(12, baseCounts));
        ASSERT_EQ(0, Pileup::variantAllele(15, baseCounts));
    }
}
