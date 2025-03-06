#pragma once

template< typename T, typename INDEX_TYPE, int N >
struct symmetricMatrix
{
  static constexpr INDEX_TYPE size() { return ( N*(N+1) ) / 2; }

  static inline INDEX_TYPE linearIndex( INDEX_TYPE const i, INDEX_TYPE const j ) 
  { 
    return (( i*(i+1) ) >> 1) + j; 
  }
  
  inline T & operator()( INDEX_TYPE const i, INDEX_TYPE const j ) 
  { 
    return m_data[linearIndex(i,j)];
  }
  
  inline T const & operator()( INDEX_TYPE const i, INDEX_TYPE const j ) const 
  { 
    return m_data[linearIndex(i,j)];
  }

  T m_data[ size() ];

};