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

#include <stdio.h>
#include "CArrayWrapper.hpp"

template< typename T, int N >
void print( T const (&a)[N], char const * name, int const numDigits )
{
  char format[8];
  snprintf( format, 8, "%%%d.%de", numDigits+7, numDigits );
  printf( "%-9s = { ", name );
  for( int i=0; i<N; ++i )
  {
    printf( format, a[i] );
    if( i < N-1 )
    {
      printf( ", " );
    }
  }
  printf( " }\n" );
}

template< typename T, int NROWS, int NCOLS >
void print( CArrayWrapper< T, NROWS, NCOLS > const & matrix, char const * name, int const numDigits )
{
  char format[8];
  snprintf( format, 8, "%%%d.%de", numDigits+7, numDigits );

  printf( "%-9s = \n{ \n", name );
  for( int i=0; i<NROWS; ++i )
  {
    printf( " { " );
    for( int j=0; j<NCOLS; ++j )
    {
      printf( format, matrix( i, j ) );
      if( j < NCOLS-1 )
      {
        printf( ", " );
      }
    }
    printf( " },\n" );
  }
  printf( " }\n" );
}
