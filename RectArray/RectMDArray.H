#ifndef _RECTMDARRAY_H_
#define _RECTMDARRAY_H_

#include "DBox.H"

#ifdef USE_CHOMBO
#include "BaseFab.H"
#endif

#include <cassert>
#include <functional>
#include <memory>

/**
 * @brief Multi-dimensional array over a rectangular domain
 *
 * Template class for a multi-dimensional array stored in a flat vector with column-major ordering. Supports
 * scalar and vector-valued data.
 *
 * @tparam T Type of data stored in the array
 * @tparam C Number of components per grid point (default 1 for scalar data)
 */
template <class T = double, unsigned int C = 1> class RectMDArray {
  public:
    /// Default constructor.
    RectMDArray();

    /**
     * @brief Constructs a RectMDArray over the Dbox a_box. Data is USUALLY initialized as zero.
     *
     * When using this constructor, it is recommended that the user initialize the data manually:
     * @code
     * Dbox B = Dbox(getZeros(),getOnes()*7);  //Domain
     * RectMDArray<double> R = RectMDArray<double>(B); //Construct array
     * R.setVal(0);  //Initialize array values as 0.
     * @endcode
     */
    RectMDArray(const DBox &a_box);

    /**
     * @brief Defines a default-constructed RectMDArray.
     *
     * Called by Constructor RectMDArray(const Dbox& a_box). Note that this can be called by an object that is
     * already been defined. If that is the case, you need to delete the memory from the previous definition
     * and allocate a new chunk of memory.
     * This also means that the data in the previous definition is lost.
     */
    void define(const DBox &a_box);

    /// Copy constructor
    RectMDArray(const RectMDArray<T, C> &a_srcArray);

    /// Move constructor
    RectMDArray(RectMDArray<T, C> &&a_srcArray) noexcept;

    /// Copy assignment operator
    RectMDArray<T, C> &operator=(const RectMDArray<T, C> &a_srcArray);

    /// Move assignment operator
    RectMDArray<T, C> &operator=(RectMDArray<T, C> &&a_srcArray) noexcept;

    /// Destructor.
    ~RectMDArray();

    /**
     * @brief Sets all values in a RectMDArray to a constant value.
     *
     * @param a_val A constant value.
     */
    void setVal(const T &a_val);

    /// Gets box over which array is defined.
    DBox getDBox() const { return m_box; };

    /**
     * @brief Copy on Intersection.
     *  Copy the part of *this's data which intersects with the domain of a_dest's domain.
     *  From this array to destination where domains overlap
     */
    void copyTo(RectMDArray<T, C> &a_dest) const;

    /**
     * @brief Indexing operator. only works for scalar RectMDArray objects
     *
     * Returns the data stored in *this corresponding to Point a_iv. Assertion error if a_iv is not in the
     * domain of *this.
     *
     * @param a_iv A Point in the domain of *this
     */
    inline T &operator[](const Point &a_iv);

    /**
     * @brief Indexing operator for const RectMDArrays. only works for scalar RectMDArray objects.
     */
    inline const T &operator[](const Point &a_iv) const;

    /**
     * @brief Indexing operator for vector-valued RectMDArray objects with a single index. Assertion failure
     * if returned type is not scalar.
     *
     * @param a_iv A Point in the domain of *this
     * @param a_comp Integer corresponding to the desired component of the data at a_iv.
     */
    inline T &operator()(const Point &a_iv, unsigned int a_comp);

    /// Indexing operator for constant vector-valued RectMDArray objects with a single index.
    inline const T &operator()(const Point &a_iv, unsigned int a_comp) const;

    /// @brief Get total number of data elements (points * components)
    inline unsigned int dataSize() { return m_box.sizeOf() * C; };

    bool defined() const;

    inline T *getPointer() const;

  private:
    DBox m_box;            ///< Domain of the array
    T *m_rawPtr = nullptr; ///< Raw data storage
};

#include "RectMDArrayImplem.H"

#endif
