#pragma once

#include "common/constants.hpp"
#include "common/CArrayWrapper.hpp"


namespace hpcReact
{
namespace bulkGeneric
{





template< typename REAL_TYPE,
          typename INT_TYPE,
          typename INDEX_TYPE,
          int NUM_SPECIES,
          int NUM_REACTIONS >
struct ParametersBase
{
  using RealType = REAL_TYPE;
  using IntType = INT_TYPE;

  static constexpr IntType numSpecies = NUM_SPECIES;
  static constexpr IntType numReactions = NUM_REACTIONS;
  RealType m_stoichiometricMatrix[numReactions][numSpecies];
};


template< typename REAL_TYPE,
          typename INT_TYPE,
          typename INDEX_TYPE,
          int NUM_REACTIONS >
struct EquilibriumKineticsModelConstants
{
  using RealType = REAL_TYPE;
  using IntType = INT_TYPE;

  static constexpr IntType numReactions = NUM_REACTIONS;

  void verifyParameters()
  {
    static constexpr int num_digits = 12;
    for( int i = 0; i < numReactions; ++i )
    {
      RealType & K = m_equilibriumConstant[i];
      RealType & kf = m_rateConstantForward[i];
      RealType & kr = m_rateConstantReverse[i];

      // Count the number of valid inputs
      int const numSpecified = (K > 0.0) + (kf > 0.0) + (kr > 0.0);

      if (numSpecified < 2)
      {
          throw std::runtime_error("Error: At least two values must be specified for reaction " + std::to_string(i));
      }
      else if (numSpecified == 2)
      {
          if (K < 0.0)       { K = kf / kr; }
          else if (kf < 0.0) { kf = K * kr; }
          else if (kr < 0.0) { kr = kf / K; }
      }
      else // numSpecified == 3
      {
        RealType const absDiff = abs( K - ( kf / kr ) );
        RealType const effectiveMagnitude = max(abs(K), abs(kf/kr));
        RealType const tolerance = effectiveMagnitude * pow(10,-num_digits);
          if ( absDiff > tolerance ) // Tolerance for floating point precision
          {
              throw std::runtime_error("Error: Inconsistent equilibrium relation for reaction " + std::to_string(i));
          }
      }
    }
  }

  RealType equilibriumConstant( int const i ) const { return m_equilibriumConstant[i]; }
  RealType rateConstantForward( int const i ) const { return m_rateConstantForward[i]; }
  RealType rateConstantReverse( int const i ) const { return m_rateConstantReverse[i]; }


  RealType m_equilibriumConstant[numReactions] = {0.0};
  RealType m_rateConstantForward[numReactions] = {0.0};
  RealType m_rateConstantReverse[numReactions] = {0.0};

};

template< typename REAL_TYPE,
          typename INT_TYPE,
          typename INDEX_TYPE,
          template < typename, typename, typename, int > typename RATE_CONSTANTS_TYPE,
          int NUM_SPECIES,
          int NUM_REACTIONS >
struct ReactionsParameters
{
  using RealType = REAL_TYPE;
  using IntType = INT_TYPE;
  using RateConstantType = RATE_CONSTANTS_TYPE< REAL_TYPE, INT_TYPE, INDEX_TYPE, NUM_REACTIONS >;
  static constexpr IntType numSpecies = NUM_SPECIES;
  static constexpr IntType numReactions = NUM_REACTIONS;
  

  static_assert( std::is_same_v< REAL_TYPE, typename RateConstantType::RealType >, "RealType of RATE_CONSTANTS_TYPE is inconsistent" );
  static_assert( std::is_same_v< INT_TYPE, typename RateConstantType::IntType >, "IntType of RATE_CONSTANTS_TYPE is inconsistent" );
  static_assert( numReactions == RateConstantType::numReactions, "numReactions of RATE_CONSTANTS_TYPE is inconsistent" );

  RealType stoichiometricMatrix( int const r, int const i ) const { return m_base.m_stoichiometricMatrix[r][i]; }
  RealType equilibriumConstant( int const r ) const { return m_rateConstants.equilibriumConstant(r); }
  RealType rateConstantForward( int const r ) const { return m_rateConstants.rateConstantForward(r); }
  RealType rateConstantReverse( int const r ) const { return m_rateConstants.rateConstantReverse(r); }

  ParametersBase< REAL_TYPE, INT_TYPE, INT_TYPE, NUM_SPECIES, NUM_REACTIONS > m_base;
  RateConstantType m_rateConstants;

  //  RealType (&m_stoichiometricMatrix)[numReactions][numSpecies] = m_base.m_stoichiometricMatrix;

};




template< typename REAL_TYPE,
          typename INT_TYPE,
          typename INDEX_TYPE,
          int NUM_PRIMARY_SPECIES,
          int NUM_KINETIC_REACTIONS >
struct KineticParameters
{
  using RealType = REAL_TYPE;
  using IntType = INT_TYPE;

  static constexpr IntType numPrimarySpecies = NUM_PRIMARY_SPECIES;
  static constexpr IntType numKineticReactions = NUM_KINETIC_REACTIONS;

  RealType m_activationEnergy[numKineticReactions];
  RealType m_equilibriumConstant[numKineticReactions];
  RealType m_stoichiometricMatrix[numKineticReactions][numPrimarySpecies];

  RealType m_rateConstant[numKineticReactions];
};

template< typename REAL_TYPE,
          typename INT_TYPE,
          typename INDEX_TYPE,
          int NUM_PRIMARY_SPECIES,
          int NUM_SECONDARY_SPECIES,
          int NUM_KINETIC_REACTIONS,
          int NUM_EQUILIBRIUM_REACTIONS >
struct BulkParameters
{
  using RealType = REAL_TYPE;
  using IntType = INT_TYPE;

  static constexpr IntType numPrimarySpecies = NUM_PRIMARY_SPECIES;
  static constexpr IntType numSecondarySpecies = NUM_SECONDARY_SPECIES;
  static constexpr IntType numKineticReactions = NUM_KINETIC_REACTIONS;
  static constexpr IntType numEquilibriumReactions = NUM_EQUILIBRIUM_REACTIONS;

};



} // namespace solidStateBattery
} // namespace hpcReact
