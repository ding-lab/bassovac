#pragma once

#include <iostream>
#include <string>

class Bassovac;
struct Sample;

class ResultFormatter {
public:
    ResultFormatter(std::ostream* out, bool fixedPoint, uint32_t precision);

    void printResult(
        const char* sequenceName,
        int32_t pos,
        int ref,
        int nVariant,
        int tVariant,
        int nBaseCounts[4],
        int tBaseCounts[4],
        const Sample& normal,
        const Sample& tumor,
        const Bassovac& bv
        );

    static std::string describeFormat();

protected:
    std::ostream* _out;
    bool _fixedPoint;
    uint32_t _precision;
};
