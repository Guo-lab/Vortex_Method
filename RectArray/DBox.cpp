#include "DBox.H"
#include <cassert>

bool reportMemory = false;
unsigned long long int memory = 0;

DBox::DBox() : m_size(0) {
    m_lowCorner = getZeros();
    m_highCorner = getOnes() * (-1); // Ensure that the box is empty
}

DBox::DBox(const DBox &a_box)
    : m_size(a_box.m_size), m_lowCorner(a_box.m_lowCorner), m_highCorner(a_box.m_highCorner) {}

void DBox::recomputeSize() {
    m_size = 1;
    for (int idir = 0; idir < DIM; idir++) {
        m_size *= size(idir);
    }
}

DBox::DBox(const Point &a_lowCorner, const Point &a_highCorner)
    : m_lowCorner(a_lowCorner), m_highCorner(a_highCorner) {
    recomputeSize();
}

// DBox returns a new box, should have put const to make clear.
DBox DBox::operator&(const DBox &a_rightDBox) const {
    Point newLow, newHigh;
    for (int i = 0; i < DIM; i++) {
        newLow[i] = std::max(m_lowCorner[i], a_rightDBox.m_lowCorner[i]);
        newHigh[i] = std::min(m_highCorner[i], a_rightDBox.m_highCorner[i]);
        // If empty intersection, return empty box
        if (newHigh[i] < newLow[i]) return DBox();
    }
    return DBox(newLow, newHigh);
}

void DBox::operator&=(const DBox &a_rightDBox) {
    DBox retval = *this & a_rightDBox;
    *this = retval;
}

DBox DBox::shift(int a_direction, int a_offset) const {
    DBox returnDBox = DBox(*this);
    returnDBox.m_lowCorner += getUnitv(a_direction) * a_offset;
    returnDBox.m_highCorner += getUnitv(a_direction) * a_offset;
    return returnDBox;
}

DBox DBox::shift(const Point &a_pt) const {
    DBox returnDBox = DBox(*this);
    returnDBox.m_lowCorner += a_pt;
    returnDBox.m_highCorner += a_pt;
    return returnDBox;
}

DBox DBox::grow(int a_offset) const {
    Point lo = m_lowCorner;
    Point hi = m_highCorner;
    lo -= getOnes() * a_offset;
    hi += getOnes() * a_offset;
    return DBox(lo, hi);
}

DBox DBox::grow(const Point &a_offset) const {
    Point lo = m_lowCorner;
    Point hi = m_highCorner;
    lo -= a_offset;
    hi += a_offset;
    return DBox(lo, hi);
}

DBox DBox::coarsen(int a_nref) const {
    assert(a_nref > 0);
    Point lo = m_lowCorner;
    Point hi = m_highCorner;
    lo /= a_nref;
    hi /= a_nref;
    return DBox(lo, hi);
}

DBox DBox::coarsen(const Point &a_pt) const {
    Point lo = m_lowCorner;
    Point hi = m_highCorner;
    lo /= a_pt;
    hi /= a_pt;
    return DBox(lo, hi);
}

DBox DBox::refine(const Point &a_pt) const {
    Point lo = m_lowCorner;
    Point hi = m_highCorner;
    lo *= a_pt;
    hi += getOnes();
    hi *= a_pt;
    hi -= getOnes();
    return DBox(lo, hi);
}
DBox DBox::refine(int a_nref) const {
    assert(a_nref > 0);
    Point lo = m_lowCorner;
    Point hi = m_highCorner;
    lo *= a_nref;
    hi += getOnes();
    hi *= a_nref;
    hi -= getOnes();
    return DBox(lo, hi);
}
DBox DBox::refineCC(const Point &a_pt) const { return refine(a_pt); }
DBox DBox::refineCC(int a_nref) const { return refine(a_nref); }

/// @brief Used in iteration.
void DBox::increment(Point &a_pt) const {
    Point current = a_pt;
    assert(DIM <= 3);

    current[0]++;

#if DIM >= 2

    if (current[0] > m_highCorner[0]) {
        current[0] = m_lowCorner[0];
        current[1]++;

#if DIM >= 3
        if (current[1] > m_highCorner[1]) {
            current[1] = m_lowCorner[1];
            current[2]++;
        }
#endif
    }

#endif

    a_pt = current;
}

Point DBox::getPoint(unsigned int k) const {
    assert(k < m_size);
    int tuple[DIM];
    for (int i = 0; i < DIM; i++) {
        int factor = (m_highCorner[i] - m_lowCorner[i] + 1);
        int kred = k % factor;
        tuple[i] = kred + m_lowCorner[i];
        k = (k - kred) / factor;
    }
    Point pt(tuple);
    return pt;
}

bool DBox::contains(const Point &a_pt) const {
    for (int idir = 0; idir < DIM; idir++) {
        if (a_pt[idir] < m_lowCorner[idir]) return false;
        if (a_pt[idir] > m_highCorner[idir]) return false;
    }
    return true;
}

void DBox::print() const { std::cout << *this << std::endl; }

Point DBox::mod(const Point &a_pt) const {
    int tuple[DIM];
    for (int i = 0; i < DIM; i++) {
        int dl = m_highCorner[i] - m_lowCorner[i] + 1;
        tuple[i] = _mod(a_pt[i], dl) + m_lowCorner[i];
    }
    return Point(tuple);
}

std::ostream &operator<<(std::ostream &os, const DBox &a_box) {
    os << "[low Corner = ";
    for (int k = 0; k < DIM; k++) {
        os << a_box.getLowCorner()[k] << " ";
    }
    os << " high Corner = ";
    for (int k = 0; k < DIM; k++) {
        os << a_box.getHighCorner()[k] << " ";
    }
    os << " size=" << a_box.sizeOf();
    os << "]";
    return os;
}

void printPoint(const Point &a_pt) {
    std::cout << "(";
    for (int dir = 0; dir < DIM; dir++) {
        std::cout << a_pt[dir];
        if (dir < DIM - 1) std::cout << ",";
    }
    std::cout << ")" << std::endl;
}

#ifdef USE_CHOMBO

DBox::DBox(const Box &a_box) {
    for (int dir = 0; dir < DIM; dir++) {
        m_lowCorner[dir] = a_box.smallEnd(dir);
        m_highCorner[dir] = a_box.bigEnd(dir);
        m_size = a_box.numPts();
    }
}

// BOX CONVERSION METHODS
Box makeBox(const DBox &a_box) {
    IntVect lo = makeIntVect(a_box.getLowCorner());
    IntVect hi = makeIntVect(a_box.getHighCorner());
    return Box(lo, hi);
}

DBox makeDBox(const Box &a_box) {
    Point lo = makePoint(a_box.smallEnd());
    Point hi = makePoint(a_box.bigEnd());
    return DBox(lo, hi);
}

#endif
