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
          typename INDEX_TYPE,
          bool LOGE_CONCENTRATION >
template< typename PARAMS_DATA,
          bool CALCULATE_DERIVATIVES,
          typename ARRAY_1D_TO_CONST,
          typename ARRAY_1D,
          typename ARRAY_2D >
HPCREACT_HOST_DEVICE inline void
KineticReactions< REAL_TYPE,
                  INT_TYPE,
                  INDEX_TYPE,
                  LOGE_CONCENTRATION
                  >::computeReactionRates_impl( RealType const &, //temperature,
                                                PARAMS_DATA const & params,
                                                ARRAY_1D_TO_CONST const & speciesConcentration,
                                                ARRAY_1D & reactionRates,
                                                ARRAY_2D & reactionRatesDerivatives )
{

  if constexpr( !CALCULATE_DERIVATIVES )
  {
    HPCREACT_UNUSED_VAR( reactionRatesDerivatives );
  }

  // loop over each reaction
  for( IntType r=0; r<PARAMS_DATA::numReactions(); ++r )
  {
    // set reaction rate to zero
    reactionRates[r] = 0.0;
    // get/calculate the forward and reverse rate constants for this reaction
    RealType const forwardRateConstant = params.rateConstantForward( r ); //* exp( -params.m_activationEnergy[r] / ( constants::R *
                                                                          // temperature ) );
    RealType const reverseRateConstant = params.rateConstantReverse( r );

    if constexpr( LOGE_CONCENTRATION )
    {
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

      reactionRates[r] = forwardRateConstant * exp( productConcForward )
                         - reverseRateConstant * exp( productConcReverse );

      if constexpr( CALCULATE_DERIVATIVES )
      {
        for( IntType i = 0; i < PARAMS_DATA::numSpecies(); ++i )
        {
          RealType const s_ri = params.stoichiometricMatrix( r, i );
          if( s_ri < 0.0 )
          {
            reactionRatesDerivatives( r, i ) = forwardRateConstant * exp( productConcForward ) * (-s_ri);
          }
          else if( s_ri > 0.0 )
          {
            reactionRatesDerivatives( r, i ) = -reverseRateConstant * exp( productConcReverse ) * s_ri;
          }
          else
          {
            reactionRatesDerivatives( r, i ) = 0.0;
          }
        }
      }
    }
    else
    {
      // variables used to build the product terms for the forward and reverse reaction rates
      RealType productConcForward = 1.0;
      RealType productConcReverse = 1.0;

      RealType dProductConcForward_dC[PARAMS_DATA::numSpecies()];
      RealType dProductConcReverse_dC[PARAMS_DATA::numSpecies()];
      for( IntType i = 0; i < PARAMS_DATA::numSpecies(); ++i )
      {
        dProductConcForward_dC[i] = 1.0;
        dProductConcReverse_dC[i] = 1.0;
      }

      // build the products for the forward and reverse reaction rates
      for( IntType i = 0; i < PARAMS_DATA::numSpecies(); ++i )
      {

        RealType const s_ri = params.stoichiometricMatrix( r, i );
        RealType const productTerm_i = speciesConcentration[i] > 1e-100 ? pow( speciesConcentration[i], fabs( s_ri ) ) : 0.0;

        if( s_ri < 0.0 )
        {
          productConcForward *= productTerm_i;
        }
        else if( s_ri > 0.0 )
        {
          productConcReverse *= productTerm_i;
        }

        if constexpr( CALCULATE_DERIVATIVES )
        {

          if( s_ri < 0.0 )
          {
            for( IntType j = 0; j < PARAMS_DATA::numSpecies(); ++j )
            {
              if( i==j )
              {
                dProductConcForward_dC[j] *= -s_ri * pow( speciesConcentration[i], -s_ri-1 );
                dProductConcReverse_dC[j] = 0.0;
              }
              else
              {
                dProductConcForward_dC[j] *= productTerm_i;
              }
            }
          }
          else if( s_ri > 0.0 )
          {
            for( IntType j = 0; j < PARAMS_DATA::numSpecies(); ++j )
            {
              if( i==j )
              {
                dProductConcReverse_dC[j] *= s_ri * pow( speciesConcentration[i], s_ri-1 );
                dProductConcForward_dC[j] = 0.0;
              }
              else
              {
                dProductConcReverse_dC[j] *= productTerm_i;
              }
            }
          }
          else
          {
            dProductConcForward_dC[i] = 0.0;
            dProductConcReverse_dC[i] = 0.0;
          }
        }
      }
      reactionRates[r] = forwardRateConstant * productConcForward - reverseRateConstant * productConcReverse;

      if constexpr( CALCULATE_DERIVATIVES )
      {
        for( IntType i = 0; i < PARAMS_DATA::numSpecies(); ++i )
        {
          reactionRatesDerivatives( r, i ) = forwardRateConstant * dProductConcForward_dC[i] - reverseRateConstant * dProductConcReverse_dC[i];
        }
      }
    } // end of if constexpr ( LOGE_CONCENTRATION )
  } // end of loop over reactions
}

template< typename REAL_TYPE,
          typename INT_TYPE,
          typename INDEX_TYPE,
          bool LOGE_CONCENTRATION >
template< typename PARAMS_DATA,
          bool CALCULATE_DERIVATIVES,
          typename ARRAY_1D_TO_CONST,
          typename ARRAY_1D_SA,
          typename ARRAY_1D,
          typename ARRAY_2D >
HPCREACT_HOST_DEVICE inline void
KineticReactions< REAL_TYPE,
                  INT_TYPE,
                  INDEX_TYPE,
                  LOGE_CONCENTRATION
                  >::computeReactionRatesQuotient_impl( RealType const &, //temperature,
                                                        PARAMS_DATA const & params,
                                                        ARRAY_1D_TO_CONST const & speciesConcentration,
                                                        ARRAY_1D_SA const & surfaceArea,
                                                        ARRAY_1D & reactionRates,
                                                        ARRAY_2D & reactionRatesDerivatives )
{
  if constexpr( !CALCULATE_DERIVATIVES )
  {
    HPCREACT_UNUSED_VAR( reactionRatesDerivatives );
  }

  // loop over each reaction
  for( IntType r=0; r<PARAMS_DATA::numReactions(); ++r )
  {
    // set reaction rate to zero
    reactionRates[r] = 0.0;

    if constexpr( CALCULATE_DERIVATIVES )
    {
      for( IntType i = 0; i < PARAMS_DATA::numSpecies(); ++i )
      {
        reactionRatesDerivatives( r, i ) = 0.0;
      }
    }

    // get/calculate the forward and reverse rate constants for this reaction
    RealType const rateConstant = params.rateConstantForward( r ); //* exp( -params.m_activationEnergy[r] / ( constants::R *
    // temperature ) );
    RealType const equilibriumConstant = params.equilibriumConstant( r );

    RealType quotient = 1.0;

    if constexpr( LOGE_CONCENTRATION )
    {
      RealType logQuotient = 0.0;
      // build the products for the forward and reverse reaction rates
      for( IntType i = 0; i < PARAMS_DATA::numSpecies(); ++i )
      {
        RealType const s_ri = params.stoichiometricMatrix( r, i );
        logQuotient += s_ri * speciesConcentration[i];
      }
      quotient = exp( logQuotient );

      if constexpr( CALCULATE_DERIVATIVES )
      {
        for( IntType i = 0; i < PARAMS_DATA::numSpecies(); ++i )
        {
          RealType const s_ri = params.stoichiometricMatrix( r, i );
          reactionRatesDerivatives( r, i ) = -rateConstant * surfaceArea[r] * s_ri * quotient / equilibriumConstant;
        }
      } // end of if constexpr ( CALCULATE_DERIVATIVES )
    } // end of if constexpr ( LOGE_CONCENTRATION )
    else
    {
      for( IntType i = 0; i < PARAMS_DATA::numSpecies(); ++i )
      {

        RealType const s_ri = params.stoichiometricMatrix( r, i );
        RealType const productTerm_i = speciesConcentration[i] > 1e-100 ? pow( speciesConcentration[i], s_ri ) : 0.0;
        quotient *= productTerm_i;
      }
      
      if constexpr( CALCULATE_DERIVATIVES )
      {
        for( IntType i = 0; i < PARAMS_DATA::numSpecies(); ++i )
        {
          RealType const s_ri = params.stoichiometricMatrix( r, i );
          if( s_ri > 0.0 || s_ri < 0.0 )
          {
            reactionRatesDerivatives( r, i ) = -rateConstant * surfaceArea[r] * s_ri * quotient / ( equilibriumConstant * speciesConcentration[i] );
          }
          else
          {
            reactionRatesDerivatives( r, i ) = 0.0;
          }
        }
      } // end of if constexpr ( CALCULATE_DERIVATIVES )
    } // end of else
    reactionRates[r] = rateConstant * surfaceArea[r] * ( 1.0 - quotient / equilibriumConstant );
  }
}

// function to  the reaction rate. Includes impact of temperature, concentration, surface area, volume fraction and porosity
template< typename REAL_TYPE,
          typename INT_TYPE,
          typename INDEX_TYPE,
          bool LOGE_CONCENTRATION >
template< typename PARAMS_DATA,
          bool CALCULATE_DERIVATIVES,
          typename ARRAY_1D_TO_CONST,
          typename ARRAY_1D,
          typename ARRAY_2D >
HPCREACT_HOST_DEVICE inline void
KineticReactions< REAL_TYPE,
                  INT_TYPE,
                  INDEX_TYPE,
                  LOGE_CONCENTRATION
                  >::computeSpeciesRates_impl( RealType const & temperature,
                                               PARAMS_DATA const & params,
                                               ARRAY_1D_TO_CONST const & speciesConcentration,
                                               ARRAY_1D & speciesRates,
                                               ARRAY_2D & speciesRatesDerivatives )
{
  RealType reactionRates[PARAMS_DATA::numReactions()] = { 0.0 };
  CArrayWrapper< double, PARAMS_DATA::numReactions(), PARAMS_DATA::numSpecies() > reactionRatesDerivatives;

  if constexpr( !CALCULATE_DERIVATIVES )
  {
    HPCREACT_UNUSED_VAR( speciesRatesDerivatives );
  }

  computeReactionRates< PARAMS_DATA >( temperature, params, speciesConcentration, reactionRates, reactionRatesDerivatives );

  for( IntType i = 0; i < PARAMS_DATA::numSpecies(); ++i )
  {
    speciesRates[i] = 0.0;
    if constexpr( CALCULATE_DERIVATIVES )
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
      if constexpr( CALCULATE_DERIVATIVES )
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
          typename INDEX_TYPE,
          bool LOGE_CONCENTRATION >
template< typename PARAMS_DATA,
          typename ARRAY_1D,
          typename ARRAY_1D_TO_CONST,
          typename ARRAY_2D >
HPCREACT_HOST_DEVICE inline void
KineticReactions< REAL_TYPE,
                  INT_TYPE,
                  INDEX_TYPE,
                  LOGE_CONCENTRATION >::timeStep( RealType const dt,
                                                  RealType const & temperature,
                                                  PARAMS_DATA const & params,
                                                  ARRAY_1D_TO_CONST const & speciesConcentration_n,
                                                  ARRAY_1D & speciesConcentration,
                                                  ARRAY_1D & speciesRates,
                                                  ARRAY_2D & speciesRatesDerivatives )
{
//  constexpr int numReactions = PARAMS_DATA::numReactions();
  constexpr int numSpecies = PARAMS_DATA::numSpecies();



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
      if constexpr( LOGE_CONCENTRATION )
      {
        nonLogC = exp( speciesConcentration[i] );
        nonLogC_n = exp( speciesConcentration_n[i] );
      }
      else
      {
        nonLogC = speciesConcentration[i];
        nonLogC_n = speciesConcentration_n[i];
      }
      residual[i] = -(nonLogC - nonLogC_n - dt * speciesRates[i]);


      for( int j = 0; j < numSpecies; ++j )
      {
        speciesRatesDerivatives( i, j ) = -dt * speciesRatesDerivatives( i, j );
      }
      if constexpr( LOGE_CONCENTRATION )
      {
        speciesRatesDerivatives( i, i ) += nonLogC;
      }
      else
      {
        speciesRatesDerivatives( i, i ) += 1.0;
      }
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


//     printf( "residual = { " );
//     for( int i = 0; i < numSpecies; ++i )
//     {
//       printf( " %g, ", residual[i] );
//     }
//     printf( "}\n" );

//     printf( "Jacobian = { \n" );
//     for( int i = 0; i < numSpecies; ++i )
//     {
//       printf( " { " );
//       for( int j = 0; j < numSpecies; ++j )
//       {
//         printf( " %g ", speciesRatesDerivatives( i, j ) );
// //        printf( " %g ", speciesRatesDerivatives( i, j ) / exp(speciesConcentration[j]) );
//         if( j < numSpecies-1 )
//         {
//           printf( ", " );
//         }
//       }
//       printf( "}, \n" );
//     }
//     printf( "}\n" );

    solveNxN_pivoted< double, numSpecies >( speciesRatesDerivatives.data, residual, deltaPrimarySpeciesConcentration );

    for( int i = 0; i < numSpecies; ++i )
    {
//      printf( "species %2d: concentration = %e, residual = %e, delta = %e \n", i, speciesConcentration[i], residual[i],
// deltaPrimarySpeciesConcentration[i] );
      speciesConcentration[i] = speciesConcentration[i] + deltaPrimarySpeciesConcentration[i];
    }

  }
}
} // namespace reactionsSystems
} // namespace hpcReact

#include "common/macrosCleanup.hpp"
