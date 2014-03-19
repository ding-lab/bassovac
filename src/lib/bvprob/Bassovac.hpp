#pragma once

#include "Sample.hpp"

#include <vector>

class Bassovac {
public:
    Bassovac(
        Sample& normal,
        Sample& tumor,
        double normalHeterozygousVariantRate,
        double normalHomozygousVariantRate,
        double tumorBackgroundMutationRate
        );

    double storeJointGenotypeProbability(
        unsigned n1, unsigned n2, unsigned t1, unsigned t2);

    double piecewisePsi(AlleleType normal, AlleleType tumor) const;

    double priorProbabilityGenotypes(
        const AlleleType normal[2],
        const AlleleType tumor[2]
        ) const;

    void bernoulliProbabilityOfReference(
        Sample& s1,
        const AlleleType a1[2],
        Sample& s2,
        const AlleleType a2[2]
        ) const;

    double probabilityObservedGivenGenotypes(
        Sample& s1,
        const AlleleType a1[2],
        Sample& s2,
        const AlleleType a2[2]
        ) const;

    double homozygousVariantProbability() const;
    double heterozygousVariantProbability() const;
    double somaticVariantProbability() const;
    double lossOfHeterozygosityProbability() const;
    double nonNotableEventProbability() const;

    double probabilityOfData() const {
        return _invProbData;
    }

protected:
    Sample& _normal;
    Sample& _tumor;
    double _normalHeterozygousVariantRate;
    double _normalHomozygousVariantRate;
    double _tumorBackgroundMutationRate;
    double _invProbData;
    double _pGenotype[2][2][2][2];
};
