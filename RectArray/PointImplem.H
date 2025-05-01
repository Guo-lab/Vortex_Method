#ifndef POINT_IMPLEM_H
#define POINT_IMPLEM_H

#pragma once
#ifndef POINT_H
#include "Point.H"
#endif

#include <cassert>

Point::Point(const int a_tuple[DIM]) {
    for (int j = 0; j < DIM; j++)
        m_tuple[j] = a_tuple[j];
}

Point::Point(const std::array<int, DIM> a_tuple) : m_tuple(a_tuple) {}

Point::Point(const Point &a_pt) { *this = a_pt; }

bool Point::operator<(const Point &a_rhs) const {
    for (int i = 0; i < DIM; i++) {
        if (m_tuple[i] < a_rhs[i]) return true;
        if (m_tuple[i] > a_rhs[i]) return false;
    }
    return false;
}

Point Point::operator+(const Point &a_rhsPoint) const {
    Point result;
    for (int i = 0; i < DIM; i++)
        result[i] = m_tuple[i] + a_rhsPoint[i];
    return result;
}

Point Point::operator-(const Point &a_rhsPoint) const {
    Point result;
    for (int i = 0; i < DIM; i++)
        result[i] = m_tuple[i] - a_rhsPoint[i];
    return result;
}

Point Point::operator*(const Point &a_rhsPoint) const {
    Point result;
    for (int i = 0; i < DIM; i++)
        result[i] = m_tuple[i] * a_rhsPoint[i];
    return result;
}

Point Point::operator/(const Point &a_rhsPoint) const {
    Point result;
    for (int i = 0; i < DIM; i++) {
        assert(a_rhsPoint[i] != 0);
        result[i] = m_tuple[i] / a_rhsPoint[i];
    }
    return result;
}

Point Point::operator*(int a_nref) const {
    Point result;
    for (int i = 0; i < DIM; i++)
        result[i] = m_tuple[i] * a_nref;
    return result;
}

Point Point::operator/(int a_nref) const {
    assert(a_nref != 0);
    Point result;
    for (int i = 0; i < DIM; i++)
        result[i] = m_tuple[i] / a_nref;
    return result;
}

void Point::operator+=(const Point &a_rhsPoint) {
    for (int j = 0; j < DIM; j++)
        m_tuple[j] = m_tuple[j] + a_rhsPoint.m_tuple[j];
}

void Point::operator-=(const Point &a_rhsPoint) {
    for (int j = 0; j < DIM; j++)
        m_tuple[j] = m_tuple[j] - a_rhsPoint.m_tuple[j];
}

void Point::operator*=(const Point &a_rhsPoint) {
    for (int j = 0; j < DIM; j++)
        m_tuple[j] = m_tuple[j] * a_rhsPoint.m_tuple[j];
}

void Point::operator/=(const Point &a_rhsPoint) {
    for (int j = 0; j < DIM; j++) {
        assert(a_rhsPoint[j] != 0);
        if (m_tuple[j] * a_rhsPoint[j] >= 0) m_tuple[j] = m_tuple[j] / a_rhsPoint[j];
        else
            m_tuple[j] = -(-m_tuple[j] / a_rhsPoint[j]);
    }
}

void Point::operator+=(int a_nref) {
    for (int j = 0; j < DIM; j++)
        m_tuple[j] = m_tuple[j] + a_nref;
}

void Point::operator-=(int a_nref) {
    for (int j = 0; j < DIM; j++)
        m_tuple[j] = m_tuple[j] - a_nref;
}

void Point::operator*=(int a_nref) {
    for (int j = 0; j < DIM; j++)
        m_tuple[j] = m_tuple[j] * a_nref;
}

void Point::operator/=(int a_nref) {
    assert(a_nref != 0);
    for (int j = 0; j < DIM; j++) {
        if (m_tuple[j] * a_nref >= 0) m_tuple[j] = m_tuple[j] / a_nref;
        else
            m_tuple[j] = -(-m_tuple[j] / a_nref);
    }
}

bool Point::operator==(const Point &a_rhsPoint) const {
    for (int j = 0; j < DIM; j++) {
        if (m_tuple[j] != a_rhsPoint[j]) return false;
    }
    return true;
}

bool Point::operator!=(const Point &a_rhsPoint) const { return !(*this == a_rhsPoint); }

inline Point getOnes() {
    int ones[DIM];
    for (int k = 0; k < DIM; k++) {
        ones[k] = 1;
    }
    return Point(ones);
}

inline Point getZeros() {
    int zeros[DIM];
    for (int k = 0; k < DIM; k++) {
        zeros[k] = 0;
    }
    return Point(zeros);
}

inline Point getBasisV(int idir) {
    Point retval = getZeros();
    retval[idir] = 1;
    return retval;
}

#ifdef USE_CHOMBO

inline Point makePoint(const IntVect &a_pt) {
    Point retVal = getZeros();
    for (int ii = 0; ii < DIM; ii++) {
        retVal[ii] = a_pt[ii];
    }
    return retVal;
}

inline IntVect makeIntVect(const Point &a_pt) {
    IntVect retVal = IntVect::Zero;
    for (int ii = 0; ii < DIM; ii++) {
        retVal[ii] = a_pt[ii];
    }
    return retVal;
}

#endif

#endif // POINT_IMPLEM_H
