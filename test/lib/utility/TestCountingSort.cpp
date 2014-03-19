#include "utility/CountingSort.hpp"

#include <gtest/gtest.h>

#include <cstdint>

TEST(TestCountingSort, unsignedChar) {
    uint8_t values[5] = {255, 9, 9, 5, 0};
    countingSort<0, 255>(values, values + 5);
    EXPECT_EQ(0, values[0]);
    EXPECT_EQ(5, values[1]);
    EXPECT_EQ(9, values[2]);
    EXPECT_EQ(9, values[3]);
    EXPECT_EQ(255, values[4]);
}

TEST(TestCountingSort, signedChar) {
    int8_t values[5] = {-128, 9, 9, 5, 127};
    countingSort<-128, 127>(values, values + 5);
    EXPECT_EQ(-128, values[0]);
    EXPECT_EQ(5, values[1]);
    EXPECT_EQ(9, values[2]);
    EXPECT_EQ(9, values[3]);
    EXPECT_EQ(127, values[4]);
}
