#ifndef _RECTMDARRAYIMPELM_H_
#define _RECTMDARRAYIMPELM_H_

#pragma once
#ifndef _RECTMDARRAY_H_
#include "RectMDArray.H"
#endif

#include <algorithm>
#include <cstdlib>
#include <functional>
#include <iostream>

// Default Constructor
template <class T, unsigned int C> RectMDArray<T, C>::RectMDArray() {}

// RectMDArray::define()
template <class T, unsigned int C> void RectMDArray<T, C>::define(const DBox &a_box) {
    if (m_rawPtr) {
        delete[] m_rawPtr;
    }
    m_rawPtr = new T[a_box.sizeOf() * C];
    m_box = a_box;
}

// Construct a RectMDArray from a DBox. Calls RectMDArray::define()
template <class T, unsigned int C> RectMDArray<T, C>::RectMDArray(const DBox &a_box) : m_rawPtr(nullptr) {
    define(a_box);
}

/// Copy constructor
template <class T, unsigned int C>
inline RectMDArray<T, C>::RectMDArray(const RectMDArray<T, C> &a_srcArray) : m_rawPtr(nullptr) {
    define(a_srcArray.m_box);
    if (m_rawPtr != nullptr && a_srcArray.m_rawPtr != nullptr) {
        std::copy(a_srcArray.m_rawPtr, a_srcArray.m_rawPtr + dataSize(), m_rawPtr);
    }
}

/// Move constructor
template <class T, unsigned int C>
inline RectMDArray<T, C>::RectMDArray(RectMDArray<T, C> &&a_srcArray) noexcept
    : m_box(a_srcArray.m_box), m_rawPtr(a_srcArray.m_rawPtr) {
    // Take ownership of source's data and reset source
    a_srcArray.m_rawPtr = nullptr;
}

/// Copy assignment
template <class T, unsigned int C>
inline RectMDArray<T, C> &RectMDArray<T, C>::operator=(const RectMDArray<T, C> &a_srcArray) {
    if (this != &a_srcArray) {
        define(a_srcArray.m_box);
        if (m_rawPtr != nullptr && a_srcArray.m_rawPtr != nullptr) {
            std::copy(a_srcArray.m_rawPtr, a_srcArray.m_rawPtr + dataSize(), m_rawPtr);
        }
    }
    return *this;
}

/// Move assignment
template <class T, unsigned int C>
inline RectMDArray<T, C> &RectMDArray<T, C>::operator=(RectMDArray<T, C> &&a_srcArray) noexcept {
    if (this != &a_srcArray) {
        if (m_rawPtr != nullptr) delete[] m_rawPtr;
        m_box = a_srcArray.m_box;
        m_rawPtr = a_srcArray.m_rawPtr;
        a_srcArray.m_rawPtr = nullptr;
    }
    return *this;
}

// Destructor
template <class T, unsigned int C> RectMDArray<T, C>::~RectMDArray() {
    if (m_rawPtr != nullptr) {
        delete[] m_rawPtr;
        m_rawPtr = nullptr;
    }
}
template <class T, unsigned int C> bool RectMDArray<T, C>::defined() const { return bool(m_rawPtr); }

template <class T, unsigned int C> T *RectMDArray<T, C>::getPointer() const { return m_rawPtr; };

/// Set all values in the array to a constant value
template <class T, unsigned int C> void RectMDArray<T, C>::setVal(const T &a_val) {
    if (m_rawPtr == nullptr) return;
    for (int k = 0; k < dataSize(); k++) m_rawPtr[k] = a_val;
}

/// Copy data to another array where domains overlap
template <class T, unsigned int C> void RectMDArray<T, C>::copyTo(RectMDArray<T, C> &a_dest) const {
    if (!defined() || !a_dest.defined()) return;

    DBox intersectBox = m_box & a_dest.m_box;
    if (intersectBox.isEmpty()) return;

    for (unsigned int comp = 0; comp < C; comp++) {
        for (Point pt = intersectBox.getLowCorner(); intersectBox.notDone(pt); intersectBox.increment(pt)) {
            a_dest(pt, comp) = (*this)(pt, comp);
        }
    }
}

template <class T, unsigned int C> T &RectMDArray<T, C>::operator[](const Point &a_iv) {
    static_assert(C == 1, "operator[] only defined for scalar RectMDArray");
    assert(defined() && "Array not defined");
    assert(m_box.contains(a_iv) && "Point out of bounds");

    int k = m_box.getIndex(a_iv);
    assert(k < m_box.sizeOf() && k >= 0);

    return m_rawPtr[k];
}

template <class T, unsigned int C> const T &RectMDArray<T, C>::operator[](const Point &a_iv) const {
    static_assert(C == 1, "operator[] only defined for scalar RectMDArray");
    assert(defined() && "Array not defined");
    assert(m_box.contains(a_iv) && "Point out of bounds");

    int k = m_box.getIndex(a_iv);
    assert(k < m_box.sizeOf() && k >= 0);

    return m_rawPtr[k];
}

template <class T, unsigned int C> T &RectMDArray<T, C>::operator()(const Point &a_iv, unsigned int a_comp) {
    assert(defined() && "Array not defined");
    assert(m_box.contains(a_iv) && "Point out of bounds");
    assert(a_comp < C && "Component out of bounds");

    int pointCount = m_box.sizeOf();
    int baseIndex = m_box.getIndex(a_iv);
    assert(baseIndex < pointCount && baseIndex >= 0);

    return m_rawPtr[baseIndex + pointCount * a_comp];
}

template <class T, unsigned int C>
const T &RectMDArray<T, C>::operator()(const Point &a_iv, unsigned int a_comp) const {
    assert(defined() && "Array not defined");
    assert(m_box.contains(a_iv) && "Point out of bounds");
    assert(a_comp < C && "Component out of bounds");

    int pointCount = m_box.sizeOf();
    int baseIndex = m_box.getIndex(a_iv);
    assert(baseIndex < pointCount && baseIndex >= 0);

    return m_rawPtr[baseIndex + pointCount * a_comp];
}

#endif // _RECTMDARRAYIMPELM_H_
