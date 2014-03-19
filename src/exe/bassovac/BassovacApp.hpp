#pragma once

#include "io/BamReaderBase.hpp"
#include "io/BamIntersector.hpp"
#include "io/BamFilter.hpp"

#include <memory>
#include <string>

class Fasta;
class Pileup;
class ResultFormatter;

class BassovacApp {
public:
    BassovacApp(int& argc, char** argv);
    ~BassovacApp();

    void run();

protected:
    void resultCb(int32_t pos, const Pileup& normal, const Pileup& tumor);

    void openBams();

protected:
    std::string _fasta;
    std::string _normalBam;
    std::string _tumorBam;
    std::string _outputFile;
    std::string _bamRegionString;
    std::unique_ptr<Fasta> _refSeq;
    std::unique_ptr<BamReaderBase> _normalReader;
    std::unique_ptr<BamReaderBase> _tumorReader;
    std::unique_ptr<ResultFormatter> _formatter;
    std::unique_ptr<BamFilter> _bamFilter;

    bool _fixedPoint;
    uint32_t _fpPrecision;
    double _normalVariantFrequency;
    double _normalPurity;
    double _tumorVariantFrequency;
    double _tumorPurity;
    double _tumorMassFraction;
    double _normalHetVariantRate;
    double _normalHomVariantRate;
    double _tumorBgMutationRate;
    double _minSomaticPvalue;
    uint32_t _minMapQual;
    uint32_t _minBaseQual;
    uint32_t _maxBins;
    uint32_t _maxDepth;
};
