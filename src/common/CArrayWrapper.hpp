#pragma once
#include "macros.hpp"


/**
 * @brief A generic wrapper for C-style arrays of varying dimensionality.
 *
 * This primary template is left undefined to catch unsupported numbers of dimensions
 * via a static_assert. Specializations exist for 1D, 2D, and 3D arrays.
 *
 * @tparam T      The type of the elements stored in the array.
 * @tparam NDIMS  The number of dimensions of the array.
 * @tparam DIMS   Parameter pack specifying the size of each dimension.
 */
template< typename T, int ... DIMS >
struct CArrayWrapper;

/**
 * @brief Specialization of CArrayWrapper for 1-dimensional arrays.
 *
 * Provides operator() and operator[] for element access and stores data contiguously in a single C-style array of size DIM0.
 *
 * @tparam T     The type of the elements stored in the array.
 * @tparam DIM0  The size of the single dimension.
 */
template< typename T, int DIM0 >
struct CArrayWrapper< T, DIM0 >
{

  /**
   * @brief Read/write access to an element by index.
   * @param dim The index (must be in range [0, DIM0)).
   * @return Reference to the element at the specified index.
   */
  HPCREACT_HOST_DEVICE inline T & operator()( int const dim ) { return data[dim]; }

  /**
   * @brief Read-only access to an element by index (const overload).
   * @param dim The index (must be in range [0, DIM0)).
   * @return Const reference to the element at the specified index.
   */
  HPCREACT_HOST_DEVICE inline T const & operator()( int const dim ) const { return data[dim]; }

  /**
   * @brief Subscript operator for read/write access.
   * @param dim The index (must be in range [0, DIM0)).
   * @return Reference to the element at the specified index.
   */
  HPCREACT_HOST_DEVICE inline T & operator[]( int const dim ) { return data[dim]; }

  /**
   * @brief Subscript operator for read-only access (const overload).
   * @param dim The index (must be in range [0, DIM0)).
   * @return Const reference to the element at the specified index.
   */
  HPCREACT_HOST_DEVICE inline T const & operator[]( int const dim ) const { return data[dim]; }

  /// The underlying 1D C-style array.
  T data[DIM0];
};

/**
 * @brief Specialization of CArrayWrapper for 2-dimensional arrays.
 *
 * Provides operator() and operator[] for element and row access. Internally,
 * it stores a 2D C-style array of size [DIM0][DIM1].
 *
 * @tparam T     The type of the elements stored in the array.
 * @tparam DIM0  The size of the first dimension.
 * @tparam DIM1  The size of the second dimension.
 */
template< typename T, int DIM0, int DIM1 >
struct CArrayWrapper< T, DIM0, DIM1 >
{
  /**
   * @brief Read/write access to an element by 2D indices.
   * @param dim0 Index in the first dimension (range [0, DIM0)).
   * @param dim1 Index in the second dimension (range [0, DIM1)).
   * @return Reference to the element at the specified 2D location.
   */
  HPCREACT_HOST_DEVICE inline T & operator()( int const dim0, int const dim1 )
  {
    return data[dim0][dim1];
  }

  /**
   * @brief Read-only access to an element by 2D indices (const overload).
   * @param dim0 Index in the first dimension (range [0, DIM0)).
   * @param dim1 Index in the second dimension (range [0, DIM1)).
   * @return Const reference to the element at the specified 2D location.
   */
  HPCREACT_HOST_DEVICE inline T const & operator()( int const dim0, int const dim1 ) const
  {
    return data[dim0][dim1];
  }

  /**
   * @brief Subscript operator returning a reference to a row in the 2D array.
   * @param dim0 The row index (range [0, DIM0)).
   * @return Reference to an array of type T[DIM1].
   *
   * This allows usage like `obj[dim0][dim1]`.
   */
  HPCREACT_HOST_DEVICE inline T       ( & operator[]( int const dim0 ))[DIM1]
  {
    return data[dim0];
  }

  /**
   * @brief Subscript operator returning a const reference to a row in the 2D array (const overload).
   * @param dim0 The row index (range [0, DIM0)).
   * @return Const reference to an array of type T[DIM1].
   */
  HPCREACT_HOST_DEVICE inline T const (&operator[]( int const dim0 ) const)[DIM1]
  {
    return data[dim0];
  }

  /// The underlying 2D C-style array of size DIM0 x DIM1.
  T data[DIM0][DIM1];
};

/**
 * @brief Specialization of CArrayWrapper for 3-dimensional arrays.
 *
 * Provides operator() and operator[] for element and row access. Internally,
 * it stores a 3D C-style array of size [DIM0][DIM1][DIM2].
 *
 * @tparam T     The type of the elements stored in the array.
 * @tparam DIM0  The size of the first dimension.
 * @tparam DIM1  The size of the second dimension.
 * @tparam DIM2  The size of the third dimension.
 */
template< typename T, int DIM0, int DIM1, int DIM2 >
struct CArrayWrapper< T, DIM0, DIM1, DIM2 >
{
  /**
   * @brief Read/write access to an element by 3D indices.
   * @param dim0 Index in the first dimension (range [0, DIM0)).
   * @param dim1 Index in the second dimension (range [0, DIM1)).
   * @param dim2 Index in the third dimension (range [0, DIM2)).
   * @return Reference to the element at the specified 3D location.
   *
   * @note Currently, this function incorrectly indexes data[dim0][dim1], missing dim2.
   *       It should be `data[dim0][dim1][dim2]`. Please correct if intended.
   */
  HPCREACT_HOST_DEVICE inline T & operator()( int const dim0, int const dim1, int const dim2 )
  {
    // NOTE: This looks like a bug in your original code. Should be data[dim0][dim1][dim2].
    return data[dim0][dim1][dim2];
  }

  /**
   * @brief Read-only access to an element by 3D indices (const overload).
   * @param dim0 Index in the first dimension (range [0, DIM0)).
   * @param dim1 Index in the second dimension (range [0, DIM1)).
   * @param dim2 Index in the third dimension (range [0, DIM2)).
   * @return Const reference to the element at the specified 3D location.
   */
  HPCREACT_HOST_DEVICE inline T const & operator()( int const dim0, int const dim1, int const dim2 ) const
  {
    // NOTE: Same potential bug as above. Should be data[dim0][dim1][dim2].
    return data[dim0][dim1][dim2];
  }

  /**
   * @brief Subscript operator returning a reference to a 2D slice in the 3D array.
   * @param dim0 The slice index (range [0, DIM0)).
   * @return Reference to an array of type T[DIM1][DIM2].
   *
   * This allows usage like `obj[dim0][dim1][dim2]`.
   */
  HPCREACT_HOST_DEVICE inline T       ( & operator[]( int const dim0 ))[DIM1][DIM2]
  {
    return data[dim0];
  }

  /**
   * @brief Subscript operator returning a const reference to a 2D slice in the 3D array (const overload).
   * @param dim0 The slice index (range [0, DIM0)).
   * @return Const reference to an array of type T[DIM1][DIM2].
   */
  HPCREACT_HOST_DEVICE inline T const (&operator[]( int const dim0 ) const)[DIM1][DIM2]
  {
    return data[dim0];
  }

  /// The underlying 3D C-style array of size DIM0 x DIM1 x DIM2.
  T data[DIM0][DIM1][DIM2];
};
