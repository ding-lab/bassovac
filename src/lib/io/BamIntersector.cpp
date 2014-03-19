#include "BamIntersector.hpp"
#include "PileupBuffer.hpp"

#include <algorithm>
#include <iostream>
#include <stdexcept>

using namespace std;


BamIntersector::BamIntersector(
        BamReaderBase& readerN,
        BamReaderBase& readerT,
        callback_t cb
        )
    : _readerN(readerN)
    , _readerT(readerT)
    , _cb(cb)
    , _tid(0)
    , _pos(0)
    , _region(readerN.region())
{
}

void BamIntersector::run() {
    _pn.push(_readerN.take());
    _pt.push(_readerT.take());

    try {
        while (!_pn.empty() && !_pt.empty()) {
            PosCompare cmp = _pn.front()->cmp(*_pt.front());
            if (cmp == BEFORE) {
                _pn.clearBefore(_pt.tid(), _pt.start());
            } else if (cmp == AFTER) {
                _pt.clearBefore(_pn.tid(), _pn.start());
            } else {
                doPileup();
            }
            if (_pn.empty()) _pn.push(_readerN.take());
            if (_pt.empty()) _pt.push(_readerT.take());
        }
    } catch (...) {
        cerr << "Error:\n";
        cerr << "Normal pileup buffer position: #" << _pn.tid() << ", pos " << _pn.start() << " -> " << _pn.end() << "\n";
        cerr << "Tumor pileup buffer position: #" << _pt.tid() << ", pos " << _pt.start() << " -> " << _pt.end() << "\n";
        throw;
    }
}

void BamIntersector::doPileup() {
    while (_pn.push(_readerN.peek())) _readerN.take();
    while (_pt.push(_readerT.peek())) _readerT.take();
    int32_t begin;
    uint32_t end;

    if (_tid != _pn.tid()) // new chromosome, reset _pos
        _pos = 0;
    _tid = _pn.tid();

    if (_region && _tid != _region->tid) {
        // FIXME: better error message
        throw std::logic_error("Region limiting error");
    }

    PosCompare cmp = _pn.cmp(_pt, &begin, &end);
    _pos = max(begin, _pos);
    end = min(end, _pn.front()->end());

    if (_region) {
        _pos = max(int(_region->beg), _pos);
        end = min(unsigned(_region->end), end);
    }

    if (cmp == BEFORE) {
        _pn.clear();
    } else if (cmp == AFTER) {
        _pt.clear();
    } else {
        while (uint32_t(_pos) < end) {
            Pileup* normal = _pn.pileup(_pos);
            Pileup* tumor  = _pt.pileup(_pos);
            if (!normal->empty() && !tumor->empty())
                _cb(_pos, *normal, *tumor);
            ++_pos;
            delete normal;
            delete tumor;
        }
        _pn.clearBefore(_pn.tid(), _pos);
        _pt.clearBefore(_pt.tid(), _pos);
    }
}
