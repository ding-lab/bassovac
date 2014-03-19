#pragma once

#include "PBin.hpp"

#include <cstdint>
#include <iostream>
#include <vector>

enum AlleleType {
    VAR,
    REF
};

struct Sample {
    Sample();

    void setValues(
        unsigned totalReads,
        unsigned supportingReads,
        double variantFrequency,
        double adjustedPurity,
        double adjustedPurityComplement,
        const uint8_t* sortedPvals,
        uint32_t nPvals,
        uint32_t nBins
        );

    double piecewisePhi(const AlleleType alleles[2]) const;

    unsigned totalReads;
    unsigned supportingReads;
    double variantFrequency;
    std::vector<PBin> readErrorBins;
    uint32_t nBins;
    double adjustedPurity;
    double adjustedPurityComplement;
};

std::ostream& operator<<(std::ostream& out, const Sample& s);
std::istream& operator>>(std::istream& in, Sample& s);
