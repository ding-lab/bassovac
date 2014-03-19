#include "BassovacApp.hpp"
#include "bvprob/Bassovac.hpp"
#include "bvprob/Fasta.hpp"
#include "bvprob/PBin.hpp"
#include "bvprob/ResultFormatter.hpp"
#include "bvprob/Sample.hpp"
#include "io/BamFilter.hpp"
#include "io/BamReader.hpp"
#include "io/Pileup.hpp"
#include "io/RegionLimitedBamReader.hpp"
#include "utility/Lut.hpp"
#include "version.h"

#include <bam.h>
#include <boost/program_options.hpp>

#include <cmath>
#include <cstdint>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <vector>

using namespace std;
using namespace std::placeholders;
namespace po = boost::program_options;

BassovacApp::BassovacApp(int& argc, char** argv)
    : _fixedPoint(false)
    , _normalVariantFrequency(0.5)
    , _tumorVariantFrequency(0.5)
    , _minBaseQual(0)
    , _maxBins(2)
    , _maxDepth(1000000)
{

    po::options_description helpOpts("Help");
    helpOpts.add_options()
        ("help,h", "this message")
        ("help-format", "describe output format")
        ("version,v", "display version information");

    po::options_description requiredOpts("Required Arguments");
    requiredOpts.add_options()
        ("fasta,f", po::value<string>(&_fasta), "fasta of reference sequence")
        ("normal-bam,n", po::value<string>(&_normalBam), "sorted .bam/.sam file containing normal reads")
        ("tumor-bam,t", po::value<string>(&_tumorBam), "sorted .bam/.sam file containing tumor reads")
        ("normal-purity", po::value<double>(&_normalPurity), "normal purity")
        ("tumor-purity", po::value<double>(&_tumorPurity), "tumor purity")
        ("tumor-mass-fraction,u", po::value<double>(&_tumorMassFraction)->default_value(1.0), "tumor mass fraction")
    ;

    po::options_description optionalOpts("Optional Arguments");
    optionalOpts.add_options()
        ("region,R",
            po::value<string>(&_bamRegionString),
            "Region to call variants in (e.g., 20:15000000-20000000)")

        ("output-file,o", po::value<string>(&_outputFile), "output file (empty or - means stdout, which is the default)")
        ("bins,b", po::value<uint32_t>(&_maxBins)->default_value(2), "maximum number of p-value bins to use")
        ("min-mapqual,q", po::value<uint32_t>(&_minMapQual)->default_value(0), "minimum mapping quality for reads")
        ("min-basequal,Q", po::value<uint32_t>(&_minBaseQual)->default_value(0), "minimum base quality for bases to be considered")
        ("min-somatic-pvalue,s", po::value<double>(&_minSomaticPvalue)->default_value(0.0), "minimum somatic pvalue for output to be displayed")
//        ("normal-var-freq", po::value<double>(&_normalVariantFrequency)->default_value(0.5), "normal variant frequency")
//        ("tumor-var-freq", po::value<double>(&_tumorVariantFrequency)->default_value(0.5), "tumor variant frequency")
        ("normal-het-rate", po::value<double>(&_normalHetVariantRate)->default_value(0.001), "normal heterozygous variant rate")
        ("normal-hom-rate", po::value<double>(&_normalHomVariantRate)->default_value(0.0005, "0.0005"), "normal homozygous variant rate")
        ("tumor-bg-rate", po::value<double>(&_tumorBgMutationRate)->default_value(0.000002, "2e-6"), "tumor background mutation rate")
        ("precision,p", po::value<uint32_t>(&_fpPrecision)->default_value(6), "floating point precision of output")
        ("fixed,x", "use fixed point notation (default=scientific)")
        ("max-depth,m", po::value<uint32_t>(&_maxDepth), "maximum expected read depth at any given position (used for optimization)")
    ;

    po::options_description allOpts("All Options");
    allOpts.add(helpOpts).add(requiredOpts).add(optionalOpts);

    po::positional_options_description posOpts;
    posOpts.add("fasta", 1);
    posOpts.add("normal-bam", 1);
    posOpts.add("tumor-bam", 1);

    po::variables_map vm;
    po::store(
        po::command_line_parser(argc, argv)
            .options(allOpts)
            .positional(posOpts).run(),
        vm
    );
    po::notify(vm);

    if (vm.count("version")) {
        stringstream ss;
        ss << "bassovac version " << __g_prog_version << ", (commit " << __g_commit_hash << ")" << endl;
        throw runtime_error(ss.str());
    }

    if (vm.count("help-format"))
        throw runtime_error(ResultFormatter::describeFormat());

    if (vm.count("help")) {
        stringstream ss;
        ss << allOpts;
        throw runtime_error(ss.str());
    }

    if (vm.count("fixed"))
        _fixedPoint = true;

    vector<string> requiredArguments = { "fasta", "normal-bam", "tumor-bam", "normal-purity", "tumor-purity" };

    for (auto iter = requiredArguments.begin(); iter != requiredArguments.end(); ++iter) {
        if (!vm.count(*iter)) {
            stringstream ss;
            ss << "Error: no value for required argument: " << *iter << "!" << endl << endl << allOpts;
            throw runtime_error(ss.str());
        }
    }
}

BassovacApp::~BassovacApp() {
}

void BassovacApp::resultCb(int32_t pos, const Pileup& normal, const Pileup& tumor) {
    const char* sequenceName = _normalReader->targetName(normal[0].tid);
    int ref;
    try {
        ref = bam_nt16_table[int(_refSeq->sequence(sequenceName, pos+1))];
    } catch (const length_error& e) {
        cerr << "Pileup error, probably due to alignments hanging off the end of a sequence:\n\t" << e.what() << "\n";
        return;
    }

    uint32_t nSupporting = normal.readsMatching(ref, _minBaseQual);
    uint32_t tSupporting = tumor.readsMatching(ref, _minBaseQual);
    if (nSupporting == normal.size() && tSupporting == tumor.size()) {
        return;
    }

    // compute base counts and variant alleles
    int nBaseCounts[4];
    int tBaseCounts[4];
    normal.baseCounts(nBaseCounts);
    tumor.baseCounts(tBaseCounts);
    int nVariant = Pileup::variantAllele(ref, nBaseCounts);
    int tVariant = Pileup::variantAllele(ref, tBaseCounts);


    // TODO: have the Sample objects create the bins to avoid extra copying
    vector<uint8_t> nErr = normal.baseQualities(_minBaseQual);
    vector<uint8_t> tErr = tumor.baseQualities(_minBaseQual);
    if (nErr.empty() || tErr.empty())
        return;

    Sample nSample;
    Sample tSample;
    nSample.setValues(
        nErr.size(),
        nSupporting,
        _normalVariantFrequency,
        _normalPurity,
        _tumorMassFraction*(1.0-_normalPurity), // adjusted normal purity complement
        &nErr[0],
        nErr.size(),
        _maxBins
        );

    tSample.setValues(
        tErr.size(),
        tSupporting,
        _tumorVariantFrequency,
        _tumorMassFraction*_tumorPurity, // adjusted tumor purity
        1.0-_tumorPurity, // adjusted tumor purity complement
        &tErr[0],
        tErr.size(),
        _maxBins
        );

    Bassovac bv(nSample, tSample, _normalHetVariantRate, _normalHomVariantRate, _tumorBgMutationRate);

    if (bv.somaticVariantProbability() < _minSomaticPvalue) {
        return;
    }

    _formatter->printResult(
        sequenceName,
        pos,
        ref,
        nVariant,
        tVariant,
        nBaseCounts,
        tBaseCounts,
        nSample,
        tSample,
        bv
        );
}

void BassovacApp::openBams() {
    _bamFilter.reset(new BamFilter(BAM_DEF_MASK, _minMapQual));

    if (_bamRegionString.empty()) {
        _normalReader.reset(new BamReader(_normalBam));
        _tumorReader.reset(new BamReader(_tumorBam));
    }
    else {
        _normalReader.reset(new RegionLimitedBamReader(_normalBam, _bamRegionString.c_str()));
        _tumorReader.reset(new RegionLimitedBamReader(_tumorBam, _bamRegionString.c_str()));
    }

    _normalReader->setFilter(_bamFilter.get());
    _tumorReader->setFilter(_bamFilter.get());
    _refSeq.reset(new Fasta(_fasta));
}

void BassovacApp::run() {
    Lut::init(_maxDepth);

    openBams();

    std::ostream* out(NULL);
    if (!_outputFile.empty() && _outputFile != "-") {
        out = new ofstream(_outputFile.c_str());
        if (!*out)
            throw runtime_error("Failed to open output file " + _outputFile);
    } else {
        out = &cout;
    }

    _formatter.reset(new ResultFormatter(out, _fixedPoint, _fpPrecision));

    BamIntersector intersector(*_normalReader, *_tumorReader,
        bind(&BassovacApp::resultCb, this, _1, _2, _3));

    clock_t start(clock());
    intersector.run();
    cerr << "Main loop: " << ((clock()-start)/double(CLOCKS_PER_SEC)) << "s CPU time\n";

    if (!_outputFile.empty())
        delete out;
}

