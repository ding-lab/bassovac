#include "bvprob/Bassovac.hpp"
#include "bvprob/Sample.hpp"

#include <gtest/gtest.h>

using namespace std;

const double normHetVarRate = 0.001;
const double normHomVarRate = 0.002;
const double tumBgMutRate   = 0.003;

class TestBassovac : public testing::Test {
public:
    TestBassovac()
        : normal(10, 5, 0.001, 0.01, 0.97)
        , tumor(12, 7, 0.002, 0.02, 0.96)
        , bv(normal, tumor, normHetVarRate, normHomVarRate, tumBgMutRate)
    {
    }

    Sample normal;
    Sample tumor;
    Bassovac bv;
};

TEST_F(TestBassovac, psi) {
    ASSERT_DOUBLE_EQ(1-tumBgMutRate, bv.piecewisePsi(REF, REF));
    ASSERT_DOUBLE_EQ(tumBgMutRate, bv.piecewisePsi(REF, VAR));
    ASSERT_DOUBLE_EQ(tumBgMutRate/3.0, bv.piecewisePsi(VAR, REF));
    ASSERT_DOUBLE_EQ(1-tumBgMutRate/3.0, bv.piecewisePsi(VAR, VAR));
}

TEST_F(TestBassovac, priorProbabilityGenotypesRange) {
    for (unsigned i = 0; i < 16; ++i) {
        AlleleType na[2] = { AlleleType((i>>1)&1 ), AlleleType(i&1) }; 
        AlleleType ta[2] = { AlleleType((i>>3)&1 ), AlleleType((i>>2)&1) }; 
        double p = bv.priorProbabilityGenotypes(na, ta); 
        ASSERT_LE(0.0, p);
        ASSERT_GE(1.0, p);
    }
}

TEST_F(TestBassovac, knownResult) {
    const uint32_t normalTotal = 1249;
    const uint32_t normalSupporting = 1213;
    const double normalVariantFrequency = 0.5;
    const double normalError = 0.000927351891576106;
    const double normalPurity = 1.0;

    const uint32_t tumorTotal = 577;
    const uint32_t tumorSupporting = 560;
    const double tumorVariantFrequency = 0.5;
    const double tumorError = 0.000912051276638698;
    const double tumorPurity = 0.76;

    const double normalHetVariantRate = 0.001;
    const double normalHomVariantRate = 0.0005;
    const double tumorBgMutationRate = 2.0e-6;

    // as computed by original perl script
    const double expected = 7.9680000394369197e-78;

    Sample normal(normalTotal, normalSupporting, normalVariantFrequency, normalError, normalPurity);
    Sample tumor(tumorTotal, tumorSupporting, tumorVariantFrequency, tumorError, tumorPurity);

    Bassovac bv(normal, tumor, normalHetVariantRate, normalHomVariantRate, tumorBgMutationRate);

    ASSERT_DOUBLE_EQ(bv.somaticVariantProbability(), expected);
}
