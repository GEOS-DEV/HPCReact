/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: (BSD-3-Clause)
 *
 * Copyright (c) 2025- Lawrence Livermore National Security LLC
 * All rights reserved
 *
 * See top level LICENSE files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

#pragma once
#include "macros.hpp"
#include <initializer_list>


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
  // default constructor
  HPCREACT_HOST_DEVICE
  constexpr CArrayWrapper() = default;

  /**
   * @brief Construct a CArrayWrapper from an initializer list.
   *
   * Allows brace-initialization with a list of values:
   * @code
   * CArrayWrapper< double, 3 > arr = {1.0, 2.0, 3.0};
   * @endcode
   *
   * @param init An initializer list with exactly DIM0 elements.
   *
   * @note No runtime bounds checking is performed on the initializer size.
   */
  HPCREACT_HOST_DEVICE
  constexpr CArrayWrapper( std::initializer_list< T > init )
  {
    // static_assert(init.size() == DIM0, "Size mismatch"); // needs c++20
    int i = 0;
    for( auto const & val : init )
    {
      data[i++] = val;
    }
  }
  /**
   * @brief Copy constructor.
   * @param src The source CArrayWrapper to copy from.
   */
  HPCREACT_HOST_DEVICE
  constexpr CArrayWrapper( CArrayWrapper const & src )
  {
    for( std::size_t i = 0; i < DIM0; i++ )
    {
      data[i] = src.data[i];
    }
  }

  /**
   * @brief Read/write access to an element by index.
   * @param dim The index (must be in range [0, DIM0)).
   * @return Reference to the element at the specified index.
   */
  HPCREACT_HOST_DEVICE
  constexpr HPCREACT_HOST_DEVICE inline T & operator()( int const dim ) { return data[dim]; }

  /**
   * @brief Read-only access to an element by index (const overload).
   * @param dim The index (must be in range [0, DIM0)).
   * @return Const reference to the element at the specified index.
   */
  HPCREACT_HOST_DEVICE
  constexpr HPCREACT_HOST_DEVICE inline T const & operator()( int const dim ) const { return data[dim]; }

  /**
   * @brief Subscript operator for read/write access.
   * @param dim The index (must be in range [0, DIM0)).
   * @return Reference to the element at the specified index.
   */
  HPCREACT_HOST_DEVICE
  constexpr HPCREACT_HOST_DEVICE inline T & operator[]( int const dim ) { return data[dim]; }

  /**
   * @brief Subscript operator for read-only access (const overload).
   * @param dim The index (must be in range [0, DIM0)).
   * @return Const reference to the element at the specified index.
   */
  HPCREACT_HOST_DEVICE
  constexpr HPCREACT_HOST_DEVICE inline T const & operator[]( int const dim ) const { return data[dim]; }

  /// The underlying 1D C-style array.
  T data[DIM0]{};
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
  // default constructor
  HPCREACT_HOST_DEVICE
  constexpr CArrayWrapper() = default;

  /**
   * @brief Copy constructor.
   * @param src The source CArrayWrapper to copy from.
   */
  HPCREACT_HOST_DEVICE
  constexpr CArrayWrapper( CArrayWrapper const & src )
  {
    for( std::size_t i = 0; i < DIM0; i++ )
    {
      for( std::size_t j = 0; j < DIM1; j++ )
        data[i][j] = src.data[i][j];
    }
  }

  /**
   * @brief Construct a 2D CArrayWrapper from nested initializer lists.
   *
   * Allows brace-initialization with a matrix-like structure:
   * @code
   * CArrayWrapper<double, 2, 3> mat = {
   *   {1.0, 2.0, 3.0},
   *   {4.0, 5.0, 6.0}
   * };
   * @endcode
   *
   * @param init A nested initializer list with exactly D0 rows and D1 elements per row.
   *
   * @note No runtime bounds checking is performed on the initializer dimensions.
   */
  HPCREACT_HOST_DEVICE
  constexpr CArrayWrapper( std::initializer_list< std::initializer_list< T > > init )
  {
    // static_assert(init.size() == DIM0, "Size mismatch"); // needs c++20
    int i = 0;
    for( auto const & row : init )
    {
      // static_assert(row.size() == DIM1, "Size mismatch"); // needs c++20
      int j = 0;
      for( auto const & val : row )
      {
        data[i][j++] = val;
      }
      ++i;
    }
  }

  /**
   * @brief Read/write access to an element by 2D indices.
   * @param dim0 Index in the first dimension (range [0, DIM0)).
   * @param dim1 Index in the second dimension (range [0, DIM1)).
   * @return Reference to the element at the specified 2D location.
   */
  HPCREACT_HOST_DEVICE
  constexpr HPCREACT_HOST_DEVICE inline T & operator()( int const dim0, int const dim1 )
  {
    return data[dim0][dim1];
  }

  /**
   * @brief Read-only access to an element by 2D indices (const overload).
   * @param dim0 Index in the first dimension (range [0, DIM0)).
   * @param dim1 Index in the second dimension (range [0, DIM1)).
   * @return Const reference to the element at the specified 2D location.
   */
  HPCREACT_HOST_DEVICE
  constexpr HPCREACT_HOST_DEVICE inline T const & operator()( int const dim0, int const dim1 ) const
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
  HPCREACT_HOST_DEVICE
  constexpr HPCREACT_HOST_DEVICE inline T       ( & operator[]( int const dim0 ))[DIM1]
  {
    return data[dim0];
  }

  /**
   * @brief Subscript operator returning a const reference to a row in the 2D array (const overload).
   * @param dim0 The row index (range [0, DIM0)).
   * @return Const reference to an array of type T[DIM1].
   */
  HPCREACT_HOST_DEVICE
  constexpr HPCREACT_HOST_DEVICE inline T const (&operator[]( int const dim0 ) const)[DIM1]
  {
    return data[dim0];
  }

  /// The underlying 2D C-style array of size DIM0 x DIM1.
  T data[DIM0][DIM1]{};
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
  // default constructor
  HPCREACT_HOST_DEVICE
  constexpr CArrayWrapper() = default;

  /**
   * @brief Construct a 3D CArrayWrapper from nested initializer lists.
   *
   * Enables tensor-like initialization using triple-nested braces:
   * @code
   * CArrayWrapper<double, 2, 2, 2> cube = {
   *   {
   *     {1.0, 2.0},
   *     {3.0, 4.0}
   *   },
   *   {
   *     {5.0, 6.0},
   *     {7.0, 8.0}
   *   }
   * };
   * @endcode
   *
   * @param init A three-level nested initializer list with D0 planes, D1 rows per plane,
   *             and D2 elements per row.
   *
   * @note This constructor does not perform size validation. Incorrect initializer sizes
   *       may lead to undefined behavior.
   */
  HPCREACT_HOST_DEVICE
  constexpr CArrayWrapper( std::initializer_list< std::initializer_list< std::initializer_list< T > > > init )
  {
    // static_assert(init.size() == DIM0, "Size mismatch"); // needs c++20
    int i = 0;
    for( auto const & plane : init )
    {
      // static_assert(plane.size() == DIM1, "Size mismatch"); // needs c++20
      int j = 0;
      for( auto const & row : plane )
      {
        // static_assert(row.size() == DIM2, "Size mismatch"); // needs c++20
        int k = 0;
        for( auto const & val : row )
        {
          data[i][j][k++] = val;
        }
        ++j;
      }
      ++i;
    }
  }

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
  HPCREACT_HOST_DEVICE
  constexpr HPCREACT_HOST_DEVICE inline T & operator()( int const dim0, int const dim1, int const dim2 )
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
  HPCREACT_HOST_DEVICE
  constexpr HPCREACT_HOST_DEVICE inline T const & operator()( int const dim0, int const dim1, int const dim2 ) const
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
  HPCREACT_HOST_DEVICE
  constexpr HPCREACT_HOST_DEVICE inline T ( & operator[]( int const dim0 ))[DIM1][DIM2]
  {
    return data[dim0];
  }

  /**
   * @brief Subscript operator returning a const reference to a 2D slice in the 3D array (const overload).
   * @param dim0 The slice index (range [0, DIM0)).
   * @return Const reference to an array of type T[DIM1][DIM2].
   */
  HPCREACT_HOST_DEVICE
  constexpr HPCREACT_HOST_DEVICE inline T const (&operator[]( int const dim0 ) const)[DIM1][DIM2]
  {
    return data[dim0];
  }

  /// The underlying 3D C-style array of size DIM0 x DIM1 x DIM2.
  T data[DIM0][DIM1][DIM2]{};
};
