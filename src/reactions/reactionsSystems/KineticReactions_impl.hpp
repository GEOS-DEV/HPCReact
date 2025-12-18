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

#include "common/logMath.hpp"
#include "common/constants.hpp"
#include "common/CArrayWrapper.hpp"
#include "common/DirectSystemSolve.hpp"

#include <math.h>
#include <string>
#include <iostream>

/** @file KineticReactions_impl.hpp
 *  @brief Header file for the KineticReactions implementation.
 *  @author HPC-REACT Team
 *  @date 2023
 */

namespace hpcReact
{
namespace reactionsSystems
{


template< typename REAL_TYPE,
          typename INT_TYPE,
          typename INDEX_TYPE >
template< typename PARAMS_DATA,
          bool CALCULATE_DERIVATIVES,
          typename ARRAY_1D_TO_CONST,
          typename ARRAY_1D,
          typename ARRAY_2D >
HPCREACT_HOST_DEVICE inline void
KineticReactions< REAL_TYPE,
                  INT_TYPE,
                  INDEX_TYPE
                  >::computeReactionRates_impl( RealType const &, //temperature,
                                                PARAMS_DATA const & params,
                                                ARRAY_1D_TO_CONST const & speciesConcentration,
                                                ARRAY_1D & reactionRates,
                                                ARRAY_2D & reactionRatesDerivatives )
{

  if constexpr ( !CALCULATE_DERIVATIVES )
  {
    HPCREACT_UNUSED_VAR( reactionRatesDerivatives );
  }

  // loop over each reaction
  for( IntType r=0; r<PARAMS_DATA::numReactions(); ++r )
  {
    // set reaction rate to zero
    reactionRates[r] = 0.0;
    // get/calculate the forward and reverse rate constants for this reaction
    RealType const forwardRateConstant = params.rateConstantForward( r ); //* logmath::exp( -params.m_activationEnergy[r] / ( constants::R *
                                                                          // temperature ) );
    RealType const reverseRateConstant = params.rateConstantReverse( r );

    RealType productConcForward = 0.0;
    RealType productConcReverse = 0.0;

    // build the products for the forward and reverse reaction rates
    for( IntType i = 0; i < PARAMS_DATA::numSpecies(); ++i )
    {

      RealType const s_ri = params.stoichiometricMatrix( r, i );

      if( s_ri < 0.0 )
      {
        productConcForward += (-s_ri) * speciesConcentration[i];
      }
      else if( s_ri > 0.0 )
      {
        productConcReverse += s_ri * speciesConcentration[i];
      }
    }

    reactionRates[r] = forwardRateConstant * logmath::exp( productConcForward )
                       - reverseRateConstant * logmath::exp( productConcReverse );

    if constexpr ( CALCULATE_DERIVATIVES )
    {
      for( IntType i = 0; i < PARAMS_DATA::numSpecies(); ++i )
      {
        RealType const s_ri = params.stoichiometricMatrix( r, i );
        if( s_ri < 0.0 )
        {
          reactionRatesDerivatives( r, i ) = forwardRateConstant * logmath::exp( productConcForward ) * (-s_ri);
        }
        else if( s_ri > 0.0 )
        {
          reactionRatesDerivatives( r, i ) = -reverseRateConstant * logmath::exp( productConcReverse ) * s_ri;
        }
        else
        {
          reactionRatesDerivatives( r, i ) = 0.0;
        }
      }
    }
  } // end of loop over reactions
}

template< typename REAL_TYPE,
          typename INT_TYPE,
          typename INDEX_TYPE >
template< typename PARAMS_DATA,
          bool CALCULATE_DERIVATIVES,
          typename ARRAY_1D_TO_CONST,
          typename ARRAY_1D_SA,
          typename ARRAY_1D,
          typename ARRAY_2D >
HPCREACT_HOST_DEVICE inline void
KineticReactions< REAL_TYPE,
                  INT_TYPE,
                  INDEX_TYPE
                  >::computeReactionRatesQuotient_impl( RealType const &, //temperature,
                                                        PARAMS_DATA const & params,
                                                        ARRAY_1D_TO_CONST const & speciesConcentration,
                                                        ARRAY_1D_SA const & surfaceArea,
                                                        ARRAY_1D & reactionRates,
                                                        ARRAY_2D & reactionRatesDerivatives )
{
  if constexpr ( !CALCULATE_DERIVATIVES )
  {
    HPCREACT_UNUSED_VAR( reactionRatesDerivatives );
  }

  // loop over each reaction
  for( IntType r=0; r<PARAMS_DATA::numReactions(); ++r )
  {
    // set reaction rate to zero
    reactionRates[r] = 0.0;

    if constexpr ( CALCULATE_DERIVATIVES )
    {
      for( IntType i = 0; i < PARAMS_DATA::numSpecies(); ++i )
      {
        reactionRatesDerivatives( r, i ) = 0.0;
      }
    }

    // get/calculate the forward and reverse rate constants for this reaction
    RealType const rateConstant = params.rateConstantForward( r ); //* logmath::exp( -params.m_activationEnergy[r] / ( constants::R *
    // temperature ) );
    RealType const equilibriumConstant = params.equilibriumConstant( r );

    RealType quotient = 1.0;

    RealType logQuotient = 0.0;
    // build the products for the forward and reverse reaction rates
    for( IntType i = 0; i < PARAMS_DATA::numSpecies(); ++i )
    {
      RealType const s_ri = params.stoichiometricMatrix( r, i );
      logQuotient += s_ri * speciesConcentration[i];
    }
    quotient = logmath::exp( logQuotient );

    if constexpr ( CALCULATE_DERIVATIVES )
    {
      for( IntType i = 0; i < PARAMS_DATA::numSpecies(); ++i )
      {
        RealType const s_ri = params.stoichiometricMatrix( r, i );
        reactionRatesDerivatives( r, i ) = -rateConstant * surfaceArea[r] * s_ri * quotient / equilibriumConstant;
      }
    } // end of if constexpr ( CALCULATE_DERIVATIVES )
    reactionRates[r] = rateConstant * surfaceArea[r] * ( 1.0 - quotient / equilibriumConstant );
  }
}

// function to  the reaction rate. Includes impact of temperature, concentration, surface area, volume fraction and porosity
template< typename REAL_TYPE,
          typename INT_TYPE,
          typename INDEX_TYPE >
template< typename PARAMS_DATA,
          bool CALCULATE_DERIVATIVES,
          typename ARRAY_1D_TO_CONST,
          typename ARRAY_1D,
          typename ARRAY_2D >
HPCREACT_HOST_DEVICE inline void
KineticReactions< REAL_TYPE,
                  INT_TYPE,
                  INDEX_TYPE
                  >::computeSpeciesRates_impl( RealType const & temperature,
                                               PARAMS_DATA const & params,
                                               ARRAY_1D_TO_CONST const & speciesConcentration,
                                               ARRAY_1D & speciesRates,
                                               ARRAY_2D & speciesRatesDerivatives )
{
  RealType reactionRates[PARAMS_DATA::numReactions()] = { 0.0 };
  CArrayWrapper< double, PARAMS_DATA::numReactions(), PARAMS_DATA::numSpecies() > reactionRatesDerivatives;

  if constexpr ( !CALCULATE_DERIVATIVES )
  {
    HPCREACT_UNUSED_VAR( speciesRatesDerivatives );
  }

  computeReactionRates< PARAMS_DATA >( temperature, params, speciesConcentration, reactionRates, reactionRatesDerivatives );

  for( IntType i = 0; i < PARAMS_DATA::numSpecies(); ++i )
  {
    speciesRates[i] = 0.0;
    if constexpr ( CALCULATE_DERIVATIVES )
    {
      for( IntType j = 0; j < PARAMS_DATA::numSpecies(); ++j )
      {
        speciesRatesDerivatives( i, j ) = 0.0;
      }
    }
    for( IntType r=0; r<PARAMS_DATA::numReactions(); ++r )
    {
      RealType const s_ir = params.stoichiometricMatrix( r, i );
      speciesRates[i] += s_ir * reactionRates[r];
      if constexpr ( CALCULATE_DERIVATIVES )
      {
        for( IntType j = 0; j < PARAMS_DATA::numSpecies(); ++j )
        {
          speciesRatesDerivatives( i, j ) += s_ir * reactionRatesDerivatives( r, j );
        }
      }
    }
  }
}

template< typename REAL_TYPE,
          typename INT_TYPE,
          typename INDEX_TYPE >
template< typename PARAMS_DATA,
          typename ARRAY_1D,
          typename ARRAY_1D_TO_CONST,
          typename ARRAY_2D >
HPCREACT_HOST_DEVICE inline void
KineticReactions< REAL_TYPE,
                  INT_TYPE,
                  INDEX_TYPE >::timeStep( RealType const dt,
                                          RealType const & temperature,
                                          PARAMS_DATA const & params,
                                          ARRAY_1D_TO_CONST const & speciesConcentration_n,
                                          ARRAY_1D & speciesConcentration,
                                          ARRAY_1D & speciesRates,
                                          ARRAY_2D & speciesRatesDerivatives )
{
//  static constexpr int numReactions = PARAMS_DATA::numReactions();
  static constexpr int numSpecies = PARAMS_DATA::numSpecies();



  REAL_TYPE residualNorm = 0.0;
  for( int k=0; k<20; ++k ) // newton loop
  {
//    printf( "iteration %2d: \n", k );

    computeSpeciesRates( temperature,
                         params,
                         speciesConcentration,
                         speciesRates,
                         speciesRatesDerivatives );

    double residual[numSpecies] = { 0.0 };
    double deltaPrimarySpeciesConcentration[numSpecies] = { 0.0 };

    // form residual and Jacobian
    for( int i = 0; i < numSpecies; ++i )
    {


      RealType nonLogC;
      RealType nonLogC_n;
      nonLogC = logmath::exp( speciesConcentration[i] );
      nonLogC_n = logmath::exp( speciesConcentration_n[i] );
      residual[i] = -(nonLogC - nonLogC_n - dt * speciesRates[i]);


      for( int j = 0; j < numSpecies; ++j )
      {
        speciesRatesDerivatives( i, j ) = -dt * speciesRatesDerivatives( i, j );
      }
      speciesRatesDerivatives( i, i ) += nonLogC;
    }



    residualNorm = 0.0;
    for( int j = 0; j < numSpecies; ++j )
    {
      residualNorm += residual[j] * residual[j];
    }
    residualNorm = sqrt( residualNorm );
    if( residualNorm < 1.0e-14 )
    {
      break;
    }

    solveNxN_pivoted< double, numSpecies >( speciesRatesDerivatives.data, residual, deltaPrimarySpeciesConcentration );

    for( int i = 0; i < numSpecies; ++i )
    {
      speciesConcentration[i] = speciesConcentration[i] + deltaPrimarySpeciesConcentration[i];
    }
  }
}
} // namespace reactionsSystems
} // namespace hpcReact

#include "common/macrosCleanup.hpp"
