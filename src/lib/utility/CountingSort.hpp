#pragma once

#include <cassert>
#include <cstddef>
#include <cstdint>
#include <iterator>
#include <limits>
#include <type_traits>

// no constexpr for gcc 4.4, which must be supported unfortunately.
// commenting the general version for now.
/*
template<
        typename IterType,
        typename ValueType = typename std::iterator_traits<IterType>::value_type,
        size_t maxValue = std::numeric_limits<ValueType>::max(),
        class Enable = typename std::enable_if<std::is_unsigned<ValueType>::value>::type
        >
void countingSort(IterType begin, IterType end) {
    uint32_t counts[maxValue + 1] = {0};
    auto out = begin;

    ValueType vmin(maxValue);
    ValueType vmax(0);

    while (begin != end) {
        vmin = std::min(vmin, *begin);
        vmax = std::max(vmax, *begin);
        ++counts[*begin];
        ++begin;
    }
    for (size_t i = vmin; i <= vmax; ++i) {
        for (size_t j = 0; j < counts[i]; ++j) {
            *out++ = i;
        }
    }
}
*/

// and implementing the specific thing needed by bassovac
template<ssize_t minValue, ssize_t maxValue, typename IterType>
void countingSort(IterType begin, IterType end) {
    typedef typename std::iterator_traits<IterType>::value_type ValueType;
    uint32_t counts[maxValue - minValue + 1] = {0};

    auto out = begin;

    ValueType vmin(maxValue);
    ValueType vmax(minValue);

    while (begin != end) {
        assert(*begin <= maxValue && *begin >= minValue);
        vmin = std::min(vmin, *begin);
        vmax = std::max(vmax, *begin);
        ++counts[*begin - minValue];
        ++begin;
    }

    for (ssize_t i = vmin; i <= vmax; ++i) {
        for (size_t j = 0; j < counts[i - minValue]; ++j) {
            *out++ = i;
        }
    }
}
