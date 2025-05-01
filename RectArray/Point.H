#ifndef POINT_H
#define POINT_H

#include <array>
#include <cassert>
#include <iostream>

#ifdef USE_CHOMBO
#include "IntVect.H"
#endif

#ifndef DIM
constexpr int DIM = 2;
#endif

/**
 * @brief A point in integer-valued DIM-dimensional space (arbitrary dimension)
 *
 * Represents a point in a rectangular grid with integer coordinates.
 * Supports common algebraic operations.
 */
class Point {
  public:
    /// Default Constructor
    inline Point() = default;

    /// Construct a Point using a C-style int array
    inline Point(const int a_tuple[DIM]);

    /// Construct a Point using a std::array of ints.
    inline Point(const std::array<int, DIM> a_tuple);

    /// Copy Constructor
    inline Point(const Point &a_pt);

    /// Move constructor
    inline Point(Point &&a_pt) noexcept = default;

    /// Assignment operator
    inline Point &operator=(const Point &a_pt) = default;
    inline Point &operator=(Point &&a_pt) noexcept = default;

    /// Definition of "<" operator on Points. Returns true if this* is "less than" *a_rhs.
    /// Only as ordering operator (for maps/sets)
    inline bool operator<(const Point &a_rhs) const;

    /// Componentwise arithmetic operations.
    /// Quotients are rounded down. Division of any component by 0 yields an error.
    inline Point operator+(const Point &a_rhsPoint) const;
    inline Point operator-(const Point &a_rhsPoint) const;
    inline Point operator*(const Point &a_pt) const;
    inline Point operator/(const Point &a_pt) const;
    inline Point operator*(int a_nref) const;
    inline Point operator/(int a_nref) const;

    /// Compound assignment operators
    inline void operator+=(const Point &a_pt);
    inline void operator-=(const Point &a_pt);
    inline void operator*=(const Point &a_pt);
    inline void operator/=(const Point &a_pt);
    inline void operator+=(int a_n);
    inline void operator-=(int a_n);
    inline void operator*=(int a_n);
    inline void operator/=(int a_n);

    /**
     * @brief Comparison operators
     * `==` Returns true if *this[\f$k\f$] == *a_pt[\f$k\f$] \f$\forall\,k\,\in\,[0,DIM)\f$
     * `!=` Returns false if *this == *a_pt.
     */
    inline bool operator==(const Point &a_pt) const;
    inline bool operator!=(const Point &a_pt) const;

    /// Get component of *this corresponding to a_index.
    inline const int &operator[](int a_index) const {
        assert(a_index < DIM && a_index >= 0);
        return m_tuple[a_index];
    };

    inline int &operator[](int a_index) {
        assert(a_index < DIM && a_index >= 0);
        return m_tuple[a_index];
    };

    /// Print out the contents of a Point. Output is formatted: (p1,p2,...)
    inline void print() const {
        std::cout << "(";
        for (int dir = 0; dir < DIM; dir++) {
            std::cout << m_tuple[dir];
            if (dir < DIM - 1) std::cout << ",";
        }
        std::cout << ")" << std::endl;
    };

    /// Print out the contents of a Point. Output is formatted: p1 p2 p3 ...
    inline void print_raw() const {
        for (int dir = 0; dir < DIM; dir++)
            std::cout << m_tuple[dir] << " ";
        std::cout << std::endl;
    };

  private:
    ///< Integer coordinates of the Point
    std::array<int, DIM> m_tuple = {};
};

/// Output stream operator for Point
inline std::ostream &operator<<(std::ostream &a_os, const Point &a_pt) {
    a_os << "(";
    for (int dir = 0; dir < DIM; dir++) {
        a_os << a_pt[dir];
        if (dir < DIM - 1) a_os << ",";
    }
    a_os << ")";
    return a_os;
}

/// Returns a Point of dimension DIM with all entries equal to 1.
inline Point getOnes();

/// Returns a Point of dimension DIM with all entries equal to 0.
inline Point getZeros();

/**
 * @brief getBasisV
 *
 * Returns a Point representing the idir-directional basis vector of dimension DIM
 * Output is a unit vector \f$\mathbf{v}: \mathbf{v}[k] = \delta_{idir}^k\f$. This function is identical to
 * getUnitv()
 */
inline Point getBasisV(int idir);

/**
 * @brief getUnitv
 * Returns a Point representing the idir-directional basis vector of dimension DIM.
 * Identical to getBasisV()
 * Output is a unit vector \f$\mathbf{v}: \mathbf{v}[k] = \delta_{idir}^k\f$. This function is identical to
 * getBasisV()
 */
inline Point getUnitv(int a_dir) {
    assert((a_dir < DIM) && (a_dir >= 0));
    int unitv[DIM];
    for (int k = 0; k < DIM; k++) {
        unitv[k] = 0;
    }
    unitv[a_dir] = 1;
    return Point(unitv);
}

#include "PointImplem.H"

#endif
