#pragma once

#include "BamReaderBase.hpp"
#include "Pileup.hpp"
#include "PileupBuffer.hpp"

#include <functional>
#include <sam.h>
#include <string>

class BamIntersector {
public:
    typedef std::function<void(int32_t, const Pileup&, const Pileup&)> callback_t;

    BamIntersector(
        BamReaderBase& readerN,
        BamReaderBase& readerT,
        callback_t cb
        );

    void run();
    void doPileup();

protected:
    BamReaderBase& _readerN;
    BamReaderBase& _readerT;
    callback_t _cb;
    int _tid;
    int32_t _pos;

    PileupBuffer _pn;
    PileupBuffer _pt;
    Region const* _region;
};
