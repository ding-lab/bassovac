#include "ResultFormatter.hpp" 
#include "bvprob/Bassovac.hpp"
#include "bvprob/Sample.hpp"

#include <bam.h>

#include <cmath>
#include <iomanip>
#include <sstream>
#include <string>

using namespace std;

ResultFormatter::ResultFormatter(std::ostream* out, bool fixedPoint, uint32_t precision)
    : _out(out)
    , _fixedPoint(fixedPoint)
    , _precision(precision)
{
}

void ResultFormatter::printResult(
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
    )
{
    *_out << 
        (_fixedPoint ? fixed : scientific) << 
        setprecision(_precision) << 
        sequenceName <<
        "\t" << pos <<
        "\t" << (pos+1) <<
        "\t" << bam_nt16_rev_table[ref] <<
        "\t" << (nVariant ? bam_nt16_rev_table[nVariant] : '.') <<
        "\t" << (tVariant ? bam_nt16_rev_table[tVariant] : '.') <<
        "\t" << nBaseCounts[0] << "," << nBaseCounts[1] << "," << nBaseCounts[2] << "," << nBaseCounts[3] <<
        "\t" << tBaseCounts[0] << "," << tBaseCounts[1] << "," << tBaseCounts[2] << "," << tBaseCounts[3] <<
        "\t" << normal.totalReads <<
        "\t" << normal.supportingReads <<
//        "\t" << normal.readErrorProbability <<
        "\t" << tumor.totalReads << 
        "\t" << tumor.supportingReads <<
//        "\t" << tumor.readErrorProbability <<
        "\t" << bv.homozygousVariantProbability() <<
        "\t" << bv.heterozygousVariantProbability() <<
        "\t" << bv.somaticVariantProbability() <<
        "\t" << bv.lossOfHeterozygosityProbability() <<
        "\t" << bv.nonNotableEventProbability() << 
        "\n";
}

string ResultFormatter::describeFormat() {
    string fields[] = {
        "Sequence name (chromosome)",
        "Start position (0-based)",
        "End position (0-based)",
        "Reference allele",
        "Normal variant allele",
        "Tumor variant allele",
        "Normal base occurrence counts (A,C,G,T)",
        "Tumor base occurrence counts (A,C,G,T)",
        "Normal read count at this position",
        "Normal reads supporting reference at this position",
//        "Normal read error probability (harmonic mean of base err p-values)",
        "Tumor read count at this position",
        "Tumor reads supporting reference at this position",
//        "Tumor read error probability (harmonic mean of base err p-values)",
        "Probability of homozygous variant",
        "Probability of heterozygous variant",
        "Probability of somatic variant",
        "Probability of loss of heterozygosity event",
        "Probability of 'uninteresting' event"
    };

    stringstream ss;
    ss << "\nBassovac output format (all fields are tab separated):\n\n";
    
    for (unsigned i = 0; i < sizeof(fields)/sizeof(fields[0]); ++i) {
        ss << "\t" << (i+1) << ") " << fields[i] << "\n";
    }
    return ss.str();
}
