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

/**
 * @file
 * @brief This file contains the definition of a symmetric matrix.
 *       It is used to store the symmetric part of a matrix in a
 *      compact way.
 */

/**
 * @brief A symmetric matrix is a matrix that is equal to its transpose.
 *        This class stores the symmetric part of a matrix in a compact way.
 * @tparam T The type of the elements of the matrix.
 * @tparam INDEX_TYPE The type of the indices of the matrix.
 * @tparam N The size of the matrix.
 * @details
 *   The matrix is stored in a 1D array, where the elements are
 *   stored in a compact way.
 *
 *   The linear index is calculated using
 *   the formula (i*(i+1))/2 + j, where i and j are the indices of
 *   the matrix. This corrosponds to filling the lower triangular
 *   part of the matrix row by row. So each row i has i+1 elements.
 *
 *   The size of the matrix is calculated using the
 *   formula (N*(N+1))/2, where N is the size of the matrix.
 */
template< typename T, typename INDEX_TYPE, int N >
struct symmetricMatrix
{
  /**
   * @brief the storage size of the matrix
   * @return the storage size of the matrix
   */
  static constexpr INDEX_TYPE size() { return ( N*(N+1) ) / 2; }

  /**
   * @brief calculates the linear index of the element at (i,j)
   * @param i the row index
   * @param j the column index
   * @return the linear index of the element at (i,j)
   */
  static inline INDEX_TYPE linearIndex( INDEX_TYPE const i, INDEX_TYPE const j )
  {
    return (( i*(i+1) ) >> 1) + j;
  }

  /**
   * @brief non-const access operator
   * @param i the row index
   * @param j the column index
   * @return reference to the element at (i,j)
   */
  inline T & operator()( INDEX_TYPE const i, INDEX_TYPE const j )
  {
    return m_data[linearIndex( i, j )];
  }

  /**
   * @brief const access operator
   * @param i the row index
   * @param j the column index
   * @return reference to the element at (i,j)
   */
  inline T const & operator()( INDEX_TYPE const i, INDEX_TYPE const j ) const
  {
    return m_data[linearIndex( i, j )];
  }

  /// Matrix is stored in a 1D c-array.
  T m_data[ size() ];

};
