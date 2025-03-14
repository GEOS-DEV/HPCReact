#pragma once

#include "common/constants.hpp"
#include "common/CArrayWrapper.hpp"

#include <stdexcept>
#include <string>
#include <utility>

namespace hpcReact
{
namespace bulkGeneric
{



template< typename REAL_TYPE,
          typename INT_TYPE,
          typename INDEX_TYPE,
          int NUM_SPECIES,
          int NUM_REACTIONS >
struct EquilibriumReactionsParameters
{
  using RealType = REAL_TYPE;
  using IntType = INT_TYPE;
  using IndexType = INDEX_TYPE;

  static constexpr IndexType numSpecies = NUM_SPECIES;
  static constexpr IndexType numReactions = NUM_REACTIONS;

  constexpr
  EquilibriumReactionsParameters( RealType const (&stoichiometricMatrix)[numReactions][numSpecies],
                                  RealType const (&equilibriumConstant)[numReactions] ):
    EquilibriumReactionsParameters( stoichiometricMatrix,
                                    equilibriumConstant,
                                    std::make_index_sequence< NUM_REACTIONS >(),
                                    std::make_index_sequence< NUM_REACTIONS *NUM_SPECIES >() )
  {}


  RealType stoichiometricMatrix( IndexType const r, int const i ) const { return m_stoichiometricMatrix[r][i]; }
  RealType equilibriumConstant( IndexType const r ) const { return m_equilibriumConstant[r]; }

  RealType m_stoichiometricMatrix[numReactions][numSpecies];
  RealType m_equilibriumConstant[numReactions];

private:
  HPCREACT_NO_MISSING_BRACES_OPEN
  template< std::size_t ... R, std::size_t ... RxS >
  constexpr
  EquilibriumReactionsParameters( RealType const (&stoichiometricMatrix)[numReactions][numSpecies],
                                  RealType const (&equilibriumConstant)[numReactions],
                                  std::index_sequence< R... >,
                                  std::index_sequence< RxS... > ):
    m_stoichiometricMatrix{ stoichiometricMatrix[RxS/numSpecies][RxS%numSpecies] ... },
    m_equilibriumConstant{ equilibriumConstant[R] ... }
  {}
  HPCREACT_NO_MISSING_BRACES_CLOSE
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

  static constexpr IndexType numSpecies = NUM_SPECIES;
  static constexpr IndexType numReactions = NUM_REACTIONS;

  KineticReactionsParameters( RealType const (&stoichiometricMatrix)[numReactions][numSpecies],
                              RealType const (&rateConstantForward)[numReactions],
                              RealType const (&rateConstantReverse)[numReactions] ):
    KineticReactionsParameters( stoichiometricMatrix,
                                rateConstantForward,
                                rateConstantReverse,
                                std::make_index_sequence< NUM_REACTIONS >(),
                                std::make_index_sequence< NUM_REACTIONS *NUM_SPECIES >() )
  {}


  RealType stoichiometricMatrix( IndexType const r, int const i ) const { return m_stoichiometricMatrix[r][i]; }
  RealType rateConstantForward( IndexType const r ) const { return m_rateConstantForward[r]; }
  RealType rateConstantReverse( IndexType const r ) const { return m_rateConstantReverse[r]; }


  RealType m_stoichiometricMatrix[numReactions][numSpecies];
  RealType m_rateConstantForward[numReactions];
  RealType m_rateConstantReverse[numReactions];

private:
  HPCREACT_NO_MISSING_BRACES(
    template< std::size_t ... R, std::size_t ... RxS >
    KineticReactionsParameters( RealType const (&stoichiometricMatrix)[numReactions][numSpecies],
                                RealType const (&rateConstantForward)[numReactions],
                                RealType const (&rateConstantReverse)[numReactions],
                                std::index_sequence< R... >,
                                std::index_sequence< RxS... > ) :
      m_stoichiometricMatrix{ stoichiometricMatrix[RxS/numSpecies][RxS%numSpecies] ... },
    m_rateConstantForward{ rateConstantForward[R] ... },
    m_rateConstantReverse{ rateConstantReverse[R] ... }
    {}
    )

};


template< typename REAL_TYPE,
          typename INT_TYPE,
          typename INDEX_TYPE,
          int NUM_SPECIES,
          int NUM_REACTIONS >
struct MixedReactionsParameters
{
  using RealType = REAL_TYPE;
  using IntType = INT_TYPE;
  using IndexType = INDEX_TYPE;
  static constexpr IndexType numSpecies = NUM_SPECIES;
  static constexpr IndexType numReactions = NUM_REACTIONS;

  constexpr
  EquilibriumReactionsParameters< RealType, IntType, IndexType, numSpecies, numReactions >
  equilibriumReactionsParameters() const
  {
    return {m_stoichiometricMatrix, m_equilibriumConstant};
  }

  constexpr
  KineticReactionsParameters< RealType, IntType, IndexType, numSpecies, numReactions >
  kineticReactionsParameters() const
  {
    return {m_stoichiometricMatrix, m_rateConstantForward, m_rateConstantReverse};
  }

  void verifyParameterConsistency()
  {
    static constexpr int num_digits = 12;
    for( int i = 0; i < numReactions; ++i )
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

  RealType m_stoichiometricMatrix[numReactions][numSpecies];
  RealType m_equilibriumConstant[numReactions];
  RealType m_rateConstantForward[numReactions];
  RealType m_rateConstantReverse[numReactions];
};



} // namespace solidStateBattery
} // namespace hpcReact
