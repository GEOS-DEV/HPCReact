
#pragma once

#include "common/constants.hpp"
#include "common/CArrayWrapper.hpp"
#include "common/macros.hpp"

#include <stdexcept>
#include <string>
#include <utility>

namespace hpcReact
{
namespace reactionsSystems
{



template< typename REAL_TYPE,
          typename INT_TYPE,
          typename INDEX_TYPE,
          int NUM_SPECIES,
          int NUM_REACTIONS,
          int NUM_SURFACE_REACTIONS = 0 >
struct EquilibriumReactionsParameters
{
  using RealType = REAL_TYPE;
  using IntType = INT_TYPE;
  using IndexType = INDEX_TYPE;

  static constexpr IndexType numSpecies() { return NUM_SPECIES; }

  static constexpr IndexType numReactions() { return NUM_REACTIONS; }

  static constexpr IndexType numSurfaceReactions() { return NUM_SURFACE_REACTIONS; }

  static constexpr IndexType numAqueousReactions() { return numReactions() - numSurfaceReactions(); }

  static constexpr IndexType numPrimarySpecies() { return numSpecies() - numReactions(); }

  static constexpr IndexType numSecondarySpecies() { return numSpecies() - numPrimarySpecies(); }


  constexpr
  EquilibriumReactionsParameters( CArrayWrapper< RealType, NUM_REACTIONS, NUM_SPECIES > const & stoichiometricMatrix,
                                  CArrayWrapper< RealType, NUM_REACTIONS > equilibriumConstant ):
    m_stoichiometricMatrix( stoichiometricMatrix ),
    m_equilibriumConstant( equilibriumConstant )
  {}


  RealType stoichiometricMatrix( IndexType const r, int const i ) const { return m_stoichiometricMatrix[r][i]; }
  RealType equilibriumConstant( IndexType const r ) const { return m_equilibriumConstant[r]; }

  CArrayWrapper< RealType, NUM_REACTIONS, NUM_SPECIES > m_stoichiometricMatrix;
  CArrayWrapper< RealType, NUM_REACTIONS > m_equilibriumConstant;
};

template< typename REAL_TYPE,
          typename INT_TYPE,
          typename INDEX_TYPE,
          int NUM_SPECIES,
          int NUM_REACTIONS >
struct KineticReactionsParameters
{
  using RealType = REAL_TYPE;
  using IntType = INT_TYPE;
  using IndexType = INDEX_TYPE;

  static constexpr IndexType numSpecies() { return NUM_SPECIES; }

  static constexpr IndexType numReactions() { return NUM_REACTIONS; }

  constexpr KineticReactionsParameters( CArrayWrapper< RealType, NUM_REACTIONS, NUM_SPECIES > const & stoichiometricMatrix,
                                        CArrayWrapper< RealType, NUM_REACTIONS > const & rateConstantForward,
                                        CArrayWrapper< RealType, NUM_REACTIONS > const & rateConstantReverse,
                                        CArrayWrapper< RealType, NUM_REACTIONS > const & equilibriumConstant ):
    m_stoichiometricMatrix( stoichiometricMatrix ),
    m_rateConstantForward( rateConstantForward ),
    m_rateConstantReverse( rateConstantReverse ),
    m_equilibiriumConstant( equilibriumConstant ) // Initialize to empty array
  {}


  RealType stoichiometricMatrix( IndexType const r, int const i ) const { return m_stoichiometricMatrix[r][i]; }
  RealType rateConstantForward( IndexType const r ) const { return m_rateConstantForward[r]; }
  RealType rateConstantReverse( IndexType const r ) const { return m_rateConstantReverse[r]; }
  RealType equilibriumConstant( IndexType const r ) const { return m_rateConstantForward[r] / m_rateConstantReverse[r]; }


  CArrayWrapper< RealType, NUM_REACTIONS, NUM_SPECIES > m_stoichiometricMatrix;
  CArrayWrapper< RealType, NUM_REACTIONS > m_rateConstantForward;
  CArrayWrapper< RealType, NUM_REACTIONS > m_rateConstantReverse;
  CArrayWrapper< RealType, NUM_REACTIONS > m_equilibiriumConstant;
};


template< typename REAL_TYPE,
          typename INT_TYPE,
          typename INDEX_TYPE,
          int NUM_SPECIES,
          int NUM_REACTIONS,
          int NUM_EQ_REACTIONS >
struct MixedReactionsParameters
{

  using RealType = REAL_TYPE;
  using IntType = INT_TYPE;
  using IndexType = INDEX_TYPE;

  constexpr MixedReactionsParameters() = default;

  constexpr MixedReactionsParameters( CArrayWrapper< RealType, NUM_REACTIONS, NUM_SPECIES > const & stoichiometricMatrix,
                                      CArrayWrapper< RealType, NUM_REACTIONS > const & equilibriumConstant,
                                      CArrayWrapper< RealType, NUM_REACTIONS > const & rateConstantForward,
                                      CArrayWrapper< RealType, NUM_REACTIONS > const & rateConstantReverse ):
    m_stoichiometricMatrix( stoichiometricMatrix ),
    m_equilibriumConstant( equilibriumConstant ),
    m_rateConstantForward( rateConstantForward ),
    m_rateConstantReverse( rateConstantReverse )
  {}

  static constexpr IndexType numReactions() { return NUM_REACTIONS; }

  static constexpr IndexType numKineticReactions() { return NUM_REACTIONS - NUM_EQ_REACTIONS; }

  static constexpr IndexType numEquilibriumReactions() { return NUM_EQ_REACTIONS; }

  static constexpr IndexType numSpecies() { return NUM_SPECIES; }

  static constexpr IndexType numPrimarySpecies() { return NUM_SPECIES - NUM_EQ_REACTIONS; }

  static constexpr IndexType numSecondarySpecies() { return NUM_EQ_REACTIONS; }

  constexpr
  EquilibriumReactionsParameters< RealType, IntType, IndexType, numSpecies(), numEquilibriumReactions() >
  equilibriumReactionsParameters() const
  {
    CArrayWrapper< RealType, numEquilibriumReactions(), numSpecies() > eqMatrix{};
    CArrayWrapper< RealType, numEquilibriumReactions() > eqConstants{};

    for( IntType i = 0; i < numEquilibriumReactions(); ++i )
    {
      for( IntType j = 0; j < numSpecies(); ++j )
      {
        eqMatrix( i, j ) = m_stoichiometricMatrix( i, j );
      }
      eqConstants( i ) = m_equilibriumConstant( i );
    }

    return { eqMatrix, eqConstants };
  }

  constexpr
  KineticReactionsParameters< RealType, IntType, IndexType, numSpecies(), numKineticReactions() >
  kineticReactionsParameters() const
  {
    CArrayWrapper< RealType, numKineticReactions(), numSpecies() > kineticMatrix{};
    CArrayWrapper< RealType, numKineticReactions() > rateConstantForward{};
    CArrayWrapper< RealType, numKineticReactions() > rateConstantReverse{};
    CArrayWrapper< RealType, numKineticReactions() > equilibriumConstant{};

    for( IndexType i = 0; i < numKineticReactions(); ++i )
    {
      for( IndexType j = 0; j < numSpecies(); ++j )
      {
        kineticMatrix( i, j ) = m_stoichiometricMatrix( numEquilibriumReactions() + i, j );
      }
      rateConstantForward( i ) = m_rateConstantForward( numEquilibriumReactions() + i );
      rateConstantReverse( i ) = m_rateConstantReverse( numEquilibriumReactions() + i );
      equilibriumConstant( i ) = m_equilibriumConstant( numEquilibriumReactions() + i );
    }

    return { kineticMatrix, rateConstantForward, rateConstantReverse, equilibriumConstant };
  }

  void verifyParameterConsistency()
  {
    static constexpr int num_digits = 12;
    for( int i = 0; i < numReactions(); ++i )
    {
      RealType & K = m_equilibriumConstant[i];
      RealType & kf = m_rateConstantForward[i];
      RealType & kr = m_rateConstantReverse[i];

      // Count the number of valid inputs
      int const numSpecified = (K > 0.0) + (kf > 0.0) + (kr > 0.0);

      if( numSpecified < 2 )
      {
        throw std::runtime_error( "Error: At least two values must be specified for reaction " + std::to_string( i ));
      }
      else if( numSpecified == 2 )
      {
        if( K < 0.0 )       { K = kf / kr; }
        else if( kf < 0.0 ) { kf = K * kr; }
        else if( kr < 0.0 ) { kr = kf / K; }
      }
      else // numSpecified == 3
      {
        RealType const absDiff = fabs( K - ( kf / kr ) );
        RealType const effectiveMagnitude = max( fabs( K ), fabs( kf/kr ));
        RealType const tolerance = effectiveMagnitude * pow( 10, -num_digits );
        if( absDiff > tolerance ) // Tolerance for floating point precision
        {
          throw std::runtime_error( "Error: Inconsistent equilibrium relation for reaction " + std::to_string( i ));
        }
      }
    }
  }

  RealType stoichiometricMatrix( IndexType const r, int const i ) const { return m_stoichiometricMatrix[r][i]; }
  RealType equilibriumConstant( IndexType const r ) const { return m_equilibriumConstant[r]; }
  RealType rateConstantForward( IndexType const r ) const { return m_rateConstantForward[r]; }
  RealType rateConstantReverse( IndexType const r ) const { return m_rateConstantReverse[r]; }

  CArrayWrapper< RealType, NUM_REACTIONS, NUM_SPECIES > m_stoichiometricMatrix;
  CArrayWrapper< RealType, NUM_REACTIONS > m_equilibriumConstant;
  CArrayWrapper< RealType, NUM_REACTIONS > m_rateConstantForward;
  CArrayWrapper< RealType, NUM_REACTIONS > m_rateConstantReverse;
};



} // namespace solidStateBattery
} // namespace hpcReact
