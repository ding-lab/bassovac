#include "PBin.hpp"
#include "utility/Lut.hpp"
#include "utility/Lut.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <map>
#include <set>
#include <stdexcept>
#include <vector>

using namespace std;

void PBin::setValues(const uint8_t* begin, const uint8_t* end)
{
    size = uint32_t(end-begin);
    harmonicMean = 0.0;
    if (size > 0) {
        while (begin != end)
            harmonicMean += Lut::phred2p_reciprocal(*begin++);
        harmonicMean = size / harmonicMean;
    }
}

uint32_t binPValues2(const uint8_t* sortedPVals, uint32_t n, std::vector<PBin>& rv) {
    std::vector<uint8_t> diffs(n - 1);
    uint8_t maxGap = 0.0;
    uint32_t maxIdx = 0;
    for (uint32_t i = 1; i < n; ++i) {
        uint8_t diff = sortedPVals[i] - sortedPVals[i - 1];
        if (diff > maxGap) {
            maxGap = diff;
            maxIdx = i;
        }
    }

    assert(rv.size() >= 2);

    if (maxGap == 0.0) {
        rv[0].setValues(sortedPVals, sortedPVals+n);
        return 1;
    }
    else {
        rv[0].setValues(sortedPVals, sortedPVals+maxIdx);
        rv[1].setValues(sortedPVals+maxIdx, sortedPVals+n);
        return 2;
    }
}

uint32_t binPValues(const uint8_t* sortedPVals, uint32_t n, uint32_t requestedBins, vector<PBin>& rv)
{
    if (requestedBins == 2) {
        return binPValues2(sortedPVals, n, rv);
    }

    uint32_t nBins = 0u;

    typedef std::map< double, std::vector<uint32_t> > MapType;
    MapType diffMap;
    for (uint32_t i = 1; i < n; ++i) {
        double diff = sortedPVals[i]-sortedPVals[i-1];
        diffMap[diff].push_back(i);
    }

    diffMap[0.0] = vector<uint32_t>(1u, 0u);

    for (auto i = diffMap.begin(); i != diffMap.end(); ++i)
        nBins += i->second.size();

    nBins = std::min(nBins, requestedBins);

    assert(nBins <= rv.size());

    if (nBins == 1) {
        rv[0].setValues(sortedPVals, sortedPVals+n);
    } else {
        set<uint32_t> binDemarcations;
        for (auto i = diffMap.rbegin(); binDemarcations.size() < nBins-1 && i != diffMap.rend(); ++i) {
            const vector<uint32_t>& v = i->second;
            for (auto j = v.begin(); binDemarcations.size() < nBins-1 && j != v.end(); ++j) {
                binDemarcations.insert(*j);
            }
        }
        uint32_t lastPos = 0u;
        uint32_t i = 0u;
        for (auto iter = binDemarcations.begin(); iter != binDemarcations.end(); ++iter) {
            rv[i++].setValues(sortedPVals+lastPos, sortedPVals+(*iter));
            lastPos = *iter;
        }
        rv[i].setValues(sortedPVals+lastPos, sortedPVals+n);
    }

    return nBins;
}
