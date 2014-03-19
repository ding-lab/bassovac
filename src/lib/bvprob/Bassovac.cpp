#include "Bassovac.hpp"
#include "PBin.hpp"
#include "utility/Binomial.hpp"

#include <boost/format.hpp>
#include <cmath>
#include <limits>
#include <stdexcept>

using boost::format;
using namespace std;

Bassovac::Bassovac(
        Sample& normal,
        Sample& tumor,
        double normalHeterozygousVariantRate,
        double normalHomozygousVariantRate,
        double tumorBackgroundMutationRate
        )
    : _normal(normal)
    , _tumor(tumor)
    , _normalHeterozygousVariantRate(normalHeterozygousVariantRate)
    , _normalHomozygousVariantRate(normalHomozygousVariantRate)
    , _tumorBackgroundMutationRate(tumorBackgroundMutationRate)
    , _invProbData(0.0)
{
    double probabilityOfData = 0.0;

    int V = int(VAR);
    int R = int(REF);

    probabilityOfData += storeJointGenotypeProbability(R, R, R, R);
    probabilityOfData += storeJointGenotypeProbability(R, R, V, V);
    probabilityOfData += storeJointGenotypeProbability(R, V, R, R);
    probabilityOfData += storeJointGenotypeProbability(V, R, R, R);
    probabilityOfData += storeJointGenotypeProbability(V, V, R, R);
    probabilityOfData += storeJointGenotypeProbability(V, V, R, V);
    probabilityOfData += storeJointGenotypeProbability(V, V, V, R);
    probabilityOfData += storeJointGenotypeProbability(V, V, V, V);

    probabilityOfData += 2 * storeJointGenotypeProbability(R, V, R, V);
    _pGenotype[V][R][V][R] = _pGenotype[R][V][R][V];

    probabilityOfData += 2 * storeJointGenotypeProbability(R, V, V, R);
    _pGenotype[V][R][R][V] = _pGenotype[R][V][V][R];

    probabilityOfData += 2 * storeJointGenotypeProbability(R, R, R, V);
    _pGenotype[R][R][V][R] = _pGenotype[R][R][R][V];

    probabilityOfData += 2 * storeJointGenotypeProbability(R, V, V, V);
    _pGenotype[V][R][V][V] = _pGenotype[R][V][V][V];

    if (probabilityOfData != 0.0) {
        _invProbData = 1.0 / probabilityOfData;
        // sometimes this number is so small that inverting it causes an overflow
        // if that is the case, just set it to the maximum value
        if (std::isinf(_invProbData))
            _invProbData = numeric_limits<double>::max();
    }
}

double Bassovac::storeJointGenotypeProbability(
        unsigned n1, unsigned n2, unsigned t1, unsigned t2)
{
    AlleleType nAlleles[2] = {AlleleType(n1), AlleleType(n2)};
    AlleleType tAlleles[2] = {AlleleType(t1), AlleleType(t2)};

    double probPrior = priorProbabilityGenotypes(nAlleles, tAlleles);
    double probNormal = probabilityObservedGivenGenotypes(
            _normal, nAlleles, _tumor, tAlleles
            );

    double probTumor = probabilityObservedGivenGenotypes(
            _tumor, tAlleles, _normal, nAlleles
            );

    double joint = probPrior * probNormal * probTumor;
    _pGenotype[n1][n2][t1][t2] = joint;
    return joint;
}



double Bassovac::piecewisePsi(AlleleType normal, AlleleType tumor) const {
    double p = _tumorBackgroundMutationRate;
    if (normal == VAR)
        p /= 3.0;
    if (normal == tumor)
        p = 1 - p;
    return p;
}

double Bassovac::priorProbabilityGenotypes(
    const AlleleType normal[2],
    const AlleleType tumor[2]
    ) const
{
    double compositePrior = 0.0;

    if (normal[0] == normal[1]) {
        if (normal[0] == REF) {
            compositePrior = 1
                - _normalHeterozygousVariantRate
                - _normalHomozygousVariantRate;
        } else {
            compositePrior = _normalHomozygousVariantRate;
        }
    } else {
        compositePrior = _normalHeterozygousVariantRate / 2.0;
    }

    compositePrior *= piecewisePsi(normal[0], tumor[0]);
    compositePrior *= piecewisePsi(normal[1], tumor[1]);

    return compositePrior;
}

void Bassovac::bernoulliProbabilityOfReference(
    Sample& s1,
    const AlleleType a1[2],
    Sample& s2,
    const AlleleType a2[2]
    ) const
{
    int nvar1 = 2 - int(a1[0]) - int(a1[1]);
    int nvar2 = 2 - int(a2[0]) - int(a2[1]);
    double pA = nvar1 * s1.adjustedPurity;
    double pB = nvar2 * s1.adjustedPurityComplement;

    vector<PBin>& bins = s1.readErrorBins;
    for (uint32_t i = 0; i < s1.nBins; ++i) {
        double err = bins[i].harmonicMean;
        bins[i].pObserveRef = 1 - (1-4.0/3.0 * err) * (pA+pB)/2.0 - err;
        if (bins[i].pObserveRef < 0 || bins[i].pObserveRef > 1) {
            throw runtime_error(str(format(
                "%1%:%2%: Probability value %3% is invalid"
                ) %__FILE__ %__LINE__ %bins[i].pObserveRef));
        }
    }
}

double Bassovac::probabilityObservedGivenGenotypes(
    Sample& s1,
    const AlleleType a1[2],
    Sample& s2,
    const AlleleType a2[2]
    ) const
{
    bernoulliProbabilityOfReference(s1, a1, s2, a2);
    vector<PBin> const& bins = s1.readErrorBins;
    double rv = 0.0;
    if (s1.nBins == 1) {
        rv = Binomial::pdf(bins[0].pObserveRef, s1.totalReads, s1.supportingReads);
    } else if (s1.nBins == 2) {
        const vector<PBin>& bins = s1.readErrorBins;
        rv = Binomial::pdfConvolve2(
            bins[0].pObserveRef, bins[1].pObserveRef,
            bins[0].size, bins[1].size,
            s1.supportingReads
            );
    } else {
        throw runtime_error("Unsupported number of quality bins");
    }

    if (rv < 0 || rv > 1) {
        throw runtime_error(str(format(
            "%1%:%2%: Probability value %3% is invalid"
            ) %__FILE__ %__LINE__ %rv));
    }

    return rv;
}

double Bassovac::homozygousVariantProbability() const {
    return _invProbData * _pGenotype[REF][REF][VAR][VAR];
}

double Bassovac::heterozygousVariantProbability() const {
    return _invProbData *
        (_pGenotype[REF][REF][VAR][REF] + _pGenotype[REF][REF][REF][VAR]);
}

double Bassovac::somaticVariantProbability() const {
    return homozygousVariantProbability() + heterozygousVariantProbability();
}

double Bassovac::lossOfHeterozygosityProbability() const {
    return _invProbData *
        (_pGenotype[REF][VAR][REF][REF] + _pGenotype[VAR][REF][REF][REF]);
}

double Bassovac::nonNotableEventProbability() const {
    return _invProbData *
        (
            _pGenotype[REF][REF][REF][REF] +
            _pGenotype[VAR][REF][VAR][REF] +
            _pGenotype[VAR][REF][REF][VAR] +
            _pGenotype[VAR][REF][VAR][VAR] +
            _pGenotype[REF][VAR][VAR][REF] +
            _pGenotype[REF][VAR][REF][VAR] +
            _pGenotype[REF][VAR][VAR][VAR] +
            _pGenotype[VAR][VAR][REF][REF] +
            _pGenotype[VAR][VAR][VAR][REF] +
            _pGenotype[VAR][VAR][REF][VAR] +
            _pGenotype[VAR][VAR][VAR][VAR]
        );
}

