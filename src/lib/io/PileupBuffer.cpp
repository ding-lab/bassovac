#include "PileupBuffer.hpp"

#include <algorithm>

using namespace std;

PileupBuffer::PileupBuffer()
    : _end(0)
{
}

PileupBuffer::~PileupBuffer() {
    clear();
}

int PileupBuffer::tid() const {
    if (empty())
        return -1;
    return _buf.front()->tid();
}

int32_t PileupBuffer::start() const {
    if (empty())
        return 0;
    return _buf.front()->start();
}

uint32_t PileupBuffer::end() const {
    return _end;
}

Pileup* PileupBuffer::pileup(uint32_t pos) const {
    Pileup* rv(Pileup::create());

    std::vector<BamEntry::PileupData> entries;
    entries.reserve(_buf.size());
    BamEntry::PileupData pileupData;
    for (auto iter = _buf.begin(); iter != _buf.end(); ++iter) {
        if (uint32_t((*iter)->start()) > pos)
            break;

        if ((*iter)->end() > pos) {
            if ((*iter)->resolveCigar(pos, pileupData)) {
                entries.push_back(pileupData);
            }
        }
    }
    rv->set(entries);

    return rv;
}

bool PileupBuffer::push(const BamEntry* b) {
    if (!b)
        return false;

    if (empty()) {
        _buf.push_back(b);
        _end = b->end();
        return true;
    }

    if (_buf.front()->cmp(*b) == OVERLAP) {
        _buf.push_back(b);
        return true;
    } else {
        return false;
    }
}

void PileupBuffer::clear() {
    for (auto i = _buf.begin(); i != _buf.end(); ++i)
        delete *i;
    _buf.clear();
    _end = 0;
}

void PileupBuffer::clearBefore(int tid, uint32_t pos) {
    auto iter = _buf.begin();
    for (; iter != _buf.end(); ++iter) {
        if ((*iter)->tid() > tid || ((*iter)->tid() == tid && (*iter)->end() > pos))
            break;
    }

    for (auto del = _buf.begin(); del != iter; ++del)
        delete *del;

    _buf.erase(_buf.begin(), iter);
    _end = 0;
    if (!_buf.empty())
        _end = _buf.front()->end();
}

PosCompare PileupBuffer::cmp(const PileupBuffer& other, int32_t* xstart, uint32_t* xend) const {
    if (tid() < other.tid())
        return BEFORE;
    if (other.tid() < tid())
        return AFTER;

    if (end() <= uint32_t(other.start()))
        return BEFORE;

    if (other.end() <= uint32_t(start()))
        return AFTER;

    *xstart = max(other.start(), start());
    *xend = min(other.end(), end());
    return OVERLAP;
}
