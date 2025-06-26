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
#include "symmetricMatrix.hpp"
#include <cmath>

namespace hpcReact
{

template< typename REAL_TYPE, int N >
bool isPositiveDefinite( REAL_TYPE const (&A)[N][N] )
{
  REAL_TYPE temp[N][N];

  // Copy A into temp to avoid modifying original
  for( int i = 0; i < N; i++ )
  {
    for( int j = 0; j < N; j++ )
    {
      temp[i][j] = A[i][j];
    }
  }

  for( int k = 0; k < N; k++ )
  {
    // Compute determinant of k-th leading principal minor using Gaussian elimination
    if( temp[k][k] <= 0 )
      return false; // Must be positive

    for( int i = k + 1; i < N; i++ )
    {
      REAL_TYPE factor = temp[i][k] / temp[k][k];
      for( int j = k; j < N; j++ )
      {
        temp[i][j] -= factor * temp[k][j];
      }
    }
  }
  return true;
}


template< typename REAL_TYPE, int N >
void solveNxN_Cholesky( REAL_TYPE const (&A)[N][N], REAL_TYPE const (&b)[N], REAL_TYPE (& x)[N] )
{
  REAL_TYPE L[N][N] = {{0}};

  // **Cholesky Decomposition**
  for( int i = 0; i < N; i++ )
  {
    for( int j = 0; j <= i; j++ )
    {
      REAL_TYPE sum = 0;
      for( int k = 0; k < j; k++ )
      {
        sum += L[i][k] * L[j][k];
      }
      if( i == j )
        L[i][j] = sqrt( A[i][i] - sum );
      else
        L[i][j] = (A[i][j] - sum) / L[j][j];
    }
  }

  // **Forward Substitution: Solve L y = b**
  //REAL_TYPE y[N];
  for( int i = 0; i < N; i++ )
  {
    x[i] = b[i];
    for( int j = 0; j < i; j++ )
      x[i] -= L[i][j] * x[j];
    x[i] /= L[i][i];
  }

  // **Backward Substitution: Solve L^T x = y**
  for( int i = N - 1; i >= 0; i-- )
  {
//    x[i] = y[i];
    for( int j = i + 1; j < N; j++ )
      x[i] -= L[j][i] * x[j];
    x[i] /= L[i][i];
  }
}



template< typename REAL_TYPE, int N >
void solveNxN_Cholesky( symmetricMatrix< REAL_TYPE, int, N > const & A,
                        REAL_TYPE const (&b)[N],
                        REAL_TYPE (& x)[N] )
{
  symmetricMatrix< REAL_TYPE, int, N > L = {0};

  // **Cholesky Decomposition**
  for( int i = 0; i < N; i++ )
  {
    for( int j = 0; j <= i; j++ )
    {
      REAL_TYPE sum = 0;
      for( int k = 0; k < j; k++ )
      {
        sum += L( i, k ) * L( j, k );
      }
      if( i == j )
      {
        L( i, j ) = sqrt( A( i, i ) - sum );
      }
      else
      {
        L( i, j ) = (A( i, j ) - sum) / L[j][j];
      }
    }
  }

  // **Forward Substitution: Solve L y = b**
  REAL_TYPE y[N];
  for( int i = 0; i < N; i++ )
  {
    y[i] = b[i];
    for( int j = 0; j < i; j++ )
    {
      y[i] -= L( i, j ) * y[j];
    }
    y[i] /= L( i, i );
  }

  // **Backward Substitution: Solve L^T x = y**
  for( int i = N - 1; i >= 0; i-- )
  {
    x[i] = y[i];
    for( int j = i + 1; j < N; j++ )
    {
      x[i] -= L( j, i ) * x[j];
    }
    x[i] /= L( i, i );
  }
}


template< typename REAL_TYPE, int N >
HPCREACT_HOST_DEVICE
void solveNxN_pivoted( REAL_TYPE (& A)[N][N], REAL_TYPE (& b)[N], REAL_TYPE (& x)[N] )
{


  int pivot[N]; // Row index tracker
  for( int i = 0; i < N; i++ )
  {
    pivot[i] = i;
  }

  // **Step 1: Forward Elimination with Pivoting**
  for( int k = 0; k < N-1; k++ )
  {
    // **Find Pivot Row**
    int max_row = k;
    REAL_TYPE max_val = fabs( A[pivot[k]][k] );
    for( int i = k + 1; i < N; i++ )
    {
      if( fabs( A[pivot[i]][k] ) > max_val )
      {
        max_val = fabs( A[pivot[i]][k] );
        max_row = i;
      }
    }

    // **Swap Rows in Pivot Array**
    if( max_row != k )
    {
      int temp = pivot[k];
      pivot[k] = pivot[max_row];
      pivot[max_row] = temp;
    }

    // **Gaussian Elimination**
    for( int i = k + 1; i < N; i++ )
    {
      REAL_TYPE factor = A[pivot[i]][k] / A[pivot[k]][k];
      for( int j = k; j < N; j++ )
      {
        A[pivot[i]][j] -= factor * A[pivot[k]][j];
      }
      b[pivot[i]] -= factor * b[pivot[k]];
    }
  }

  // **Step 2: Back-Substitution**
  for( int i = N - 1; i >= 0; --i )
  {
    x[i] = b[pivot[i]];
    for( int j = i + 1; j < N; j++ )
    {
      x[i] -= A[pivot[i]][j] * x[j];
    }
    x[i] /= A[pivot[i]][i]; // Normalize
  }

}

} // namespace hpcReact
