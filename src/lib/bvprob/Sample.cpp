#include "Sample.hpp"

#include <cstdint>

using namespace std;

Sample::Sample()
    : totalReads(0)
    , supportingReads(0)
    , variantFrequency(0.0)
    , adjustedPurity(0.0)
    , adjustedPurityComplement(0.0)
{
    readErrorBins.resize(5);
}

void Sample::setValues(
        unsigned totalReads,
        unsigned supportingReads,
        double variantFrequency,
        double adjustedPurity,
        double adjustedPurityComplement,
        const uint8_t* sortedPvals,
        uint32_t nPvals,
        uint32_t nBins
        )
{
    this->totalReads = totalReads;
    this->supportingReads = supportingReads;
    this->variantFrequency = variantFrequency;
    this->readErrorBins = readErrorBins;
    this->adjustedPurity = adjustedPurity;
    this->adjustedPurityComplement = adjustedPurityComplement;
    this->nBins = binPValues(sortedPvals, nPvals, nBins, readErrorBins);
}

double Sample::piecewisePhi(const AlleleType alleles[2]) const {
    if (alleles[0] == alleles[1])
        return alleles[0] == REF ? 0.0 : 1.0;
    else
        return variantFrequency;
}

std::ostream& operator<<(std::ostream& out, const Sample& s) {
    out << s.totalReads << "\t"
        << s.supportingReads << "\t"
        << s.variantFrequency << "\t"
        << s.readErrorBins.size() << "\t";
    for (auto iter = s.readErrorBins.begin(); iter != s.readErrorBins.end(); ++iter)
        out << iter->size << "\t" << iter->harmonicMean << "\t";
    out << s.adjustedPurity << "\t" << s.adjustedPurityComplement << "\n";
 
    return out;
}

std::istream& operator>>(std::istream& in, Sample& s) {
    size_t nBins;
    s.totalReads = 6;
    in >> s.totalReads 
        >> s.supportingReads 
        >> s.variantFrequency 
        >> nBins;
    s.readErrorBins.resize(nBins);
    for (auto iter = s.readErrorBins.begin(); iter != s.readErrorBins.end(); ++iter)
        in >> iter->size >> iter->harmonicMean;
    in >> s.adjustedPurity >> s.adjustedPurityComplement;
     
    return in;
}
