#pragma once

#include <cstdint>
#include <vector>

class PBin {
public:
    uint32_t size;
    double harmonicMean;
    double pObserveRef;

    void setValues(const uint8_t* begin, const uint8_t* end);
};

uint32_t binPValues2(const uint8_t* sortedPVals, uint32_t n, std::vector<PBin>& rv);
uint32_t binPValues(const uint8_t* sortedPVals, uint32_t n, uint32_t requestedBins, std::vector<PBin>& rv);
