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

  HPCREACT_HOST_DEVICE static constexpr IndexType numSpecies() { return NUM_SPECIES; }

  HPCREACT_HOST_DEVICE static constexpr IndexType numReactions() { return NUM_REACTIONS; }

  HPCREACT_HOST_DEVICE static constexpr IndexType numSurfaceReactions() { return NUM_SURFACE_REACTIONS; }

  HPCREACT_HOST_DEVICE static constexpr IndexType numAqueousReactions() { return numReactions() - numSurfaceReactions(); }

  HPCREACT_HOST_DEVICE static constexpr IndexType numPrimarySpecies() { return numSpecies() - numReactions(); }

  HPCREACT_HOST_DEVICE static constexpr IndexType numSecondarySpecies() { return numSpecies() - numPrimarySpecies(); }

  HPCREACT_HOST_DEVICE
  constexpr
  EquilibriumReactionsParameters( CArrayWrapper< RealType, NUM_REACTIONS, NUM_SPECIES > const & stoichiometricMatrix,
                                  CArrayWrapper< RealType, NUM_REACTIONS > equilibriumConstant,
                                  CArrayWrapper< IntType, NUM_REACTIONS > mobileSecondarySpeciesFlag ):
    m_stoichiometricMatrix( stoichiometricMatrix ),
    m_equilibriumConstant( equilibriumConstant ),
    m_mobileSecondarySpeciesFlag( mobileSecondarySpeciesFlag )
  {}


  HPCREACT_HOST_DEVICE RealType stoichiometricMatrix( IndexType const r, int const i ) const { return m_stoichiometricMatrix[r][i]; }
  HPCREACT_HOST_DEVICE RealType equilibriumConstant( IndexType const r ) const { return m_equilibriumConstant[r]; }
  HPCREACT_HOST_DEVICE IntType mobileSecondarySpeciesFlag( IndexType const r ) const { return m_mobileSecondarySpeciesFlag[r]; }

  CArrayWrapper< RealType, NUM_REACTIONS, NUM_SPECIES > m_stoichiometricMatrix;
  CArrayWrapper< RealType, NUM_REACTIONS > m_equilibriumConstant;
  CArrayWrapper< IntType, NUM_REACTIONS > m_mobileSecondarySpeciesFlag;
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

  HPCREACT_HOST_DEVICE static constexpr IndexType numSpecies() { return NUM_SPECIES; }

  HPCREACT_HOST_DEVICE static constexpr IndexType numReactions() { return NUM_REACTIONS; }

  HPCREACT_HOST_DEVICE
  constexpr KineticReactionsParameters( CArrayWrapper< RealType, NUM_REACTIONS, NUM_SPECIES > const & stoichiometricMatrix,
                                        CArrayWrapper< RealType, NUM_REACTIONS > const & rateConstantForward,
                                        CArrayWrapper< RealType, NUM_REACTIONS > const & rateConstantReverse,
                                        CArrayWrapper< RealType, NUM_REACTIONS > const & equilibriumConstant,
                                        IntType const reactionRatesUpdateOption ):
    m_stoichiometricMatrix( stoichiometricMatrix ),
    m_rateConstantForward( rateConstantForward ),
    m_rateConstantReverse( rateConstantReverse ),
    m_equilibiriumConstant( equilibriumConstant ), // Initialize to empty array
    m_reactionRatesUpdateOption( reactionRatesUpdateOption )
  {}


  HPCREACT_HOST_DEVICE RealType stoichiometricMatrix( IndexType const r, int const i ) const { return m_stoichiometricMatrix[r][i]; }
  HPCREACT_HOST_DEVICE RealType rateConstantForward( IndexType const r ) const { return m_rateConstantForward[r]; }
  HPCREACT_HOST_DEVICE RealType rateConstantReverse( IndexType const r ) const { return m_rateConstantReverse[r]; }
  HPCREACT_HOST_DEVICE RealType equilibriumConstant( IndexType const r ) const { return m_rateConstantForward[r] / m_rateConstantReverse[r]; }

  HPCREACT_HOST_DEVICE IntType reactionRatesUpdateOption() const { return m_reactionRatesUpdateOption; }

  CArrayWrapper< RealType, NUM_REACTIONS, NUM_SPECIES > m_stoichiometricMatrix;
  CArrayWrapper< RealType, NUM_REACTIONS > m_rateConstantForward;
  CArrayWrapper< RealType, NUM_REACTIONS > m_rateConstantReverse;
  CArrayWrapper< RealType, NUM_REACTIONS > m_equilibiriumConstant;

  IntType m_reactionRatesUpdateOption; // 0: forward and reverse rate. 1: quotient form.
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
                                      CArrayWrapper< RealType, NUM_REACTIONS > const & rateConstantReverse,
                                      CArrayWrapper< IntType, NUM_REACTIONS > mobileSecondarySpeciesFlag,
                                      IntType const reactionRatesUpdateOption = 1 ):
    m_stoichiometricMatrix( stoichiometricMatrix ),
    m_equilibriumConstant( equilibriumConstant ),
    m_rateConstantForward( rateConstantForward ),
    m_rateConstantReverse( rateConstantReverse ),
    m_mobileSecondarySpeciesFlag( mobileSecondarySpeciesFlag ),
    m_reactionRatesUpdateOption( reactionRatesUpdateOption )
  {}

  HPCREACT_HOST_DEVICE static constexpr IndexType numReactions() { return NUM_REACTIONS; }

  HPCREACT_HOST_DEVICE static constexpr IndexType numKineticReactions() { return NUM_REACTIONS - NUM_EQ_REACTIONS; }

  HPCREACT_HOST_DEVICE static constexpr IndexType numEquilibriumReactions() { return NUM_EQ_REACTIONS; }

  HPCREACT_HOST_DEVICE static constexpr IndexType numSpecies() { return NUM_SPECIES; }

  HPCREACT_HOST_DEVICE static constexpr IndexType numPrimarySpecies() { return NUM_SPECIES - NUM_EQ_REACTIONS; }

  HPCREACT_HOST_DEVICE static constexpr IndexType numSecondarySpecies() { return NUM_EQ_REACTIONS; }

  HPCREACT_HOST_DEVICE
  constexpr
  EquilibriumReactionsParameters< RealType, IntType, IndexType, numSpecies(), numEquilibriumReactions() >
  equilibriumReactionsParameters() const
  {
    CArrayWrapper< RealType, numEquilibriumReactions(), numSpecies() > eqMatrix{};
    CArrayWrapper< RealType, numEquilibriumReactions() > eqConstants{};
    CArrayWrapper< IntType, numEquilibriumReactions() > mobileSpeciesFlags{};

    for( IntType i = 0; i < numEquilibriumReactions(); ++i )
    {
      for( IntType j = 0; j < numSpecies(); ++j )
      {
        eqMatrix( i, j ) = m_stoichiometricMatrix( i, j );
      }
      eqConstants( i ) = m_equilibriumConstant( i );
      mobileSpeciesFlags( i ) = m_mobileSecondarySpeciesFlag( i );
    }

    return { eqMatrix, eqConstants, mobileSpeciesFlags };
  }

  HPCREACT_HOST_DEVICE
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

    return { kineticMatrix, rateConstantForward, rateConstantReverse, equilibriumConstant, m_reactionRatesUpdateOption };
  }

  HPCREACT_HOST_DEVICE
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
        if( absDiff > tolerance )  // Tolerance for floating point precision
        {
          throw std::runtime_error( "Error: Inconsistent equilibrium relation for reaction " + std::to_string( i ));
        }
      }
    }
  }

  HPCREACT_HOST_DEVICE RealType stoichiometricMatrix( IndexType const r, int const i ) const { return m_stoichiometricMatrix[r][i]; }
  HPCREACT_HOST_DEVICE RealType equilibriumConstant( IndexType const r ) const { return m_equilibriumConstant[r]; }
  HPCREACT_HOST_DEVICE RealType rateConstantForward( IndexType const r ) const { return m_rateConstantForward[r]; }
  HPCREACT_HOST_DEVICE RealType rateConstantReverse( IndexType const r ) const { return m_rateConstantReverse[r]; }

  CArrayWrapper< RealType, NUM_REACTIONS, NUM_SPECIES > m_stoichiometricMatrix;
  CArrayWrapper< RealType, NUM_REACTIONS > m_equilibriumConstant;
  CArrayWrapper< RealType, NUM_REACTIONS > m_rateConstantForward;
  CArrayWrapper< RealType, NUM_REACTIONS > m_rateConstantReverse;
  CArrayWrapper< IntType, NUM_REACTIONS > m_mobileSecondarySpeciesFlag;

  IntType m_reactionRatesUpdateOption; // 0: forward and reverse rate. 1: quotient form.
};



} // namespace solidStateBattery
} // namespace hpcReact
