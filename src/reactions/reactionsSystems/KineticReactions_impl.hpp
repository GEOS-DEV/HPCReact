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
#include "constitutive/activity/activity.hpp"

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
          typename ACTIVITY_MODEL,
          bool LOGE_CONCENTRATION >
template< typename PARAMS_DATA,
          bool CALCULATE_DERIVATIVES,
          typename ARRAY_1D_TO_CONST,
          typename ARRAY_1D,
          typename ARRAY_2D,
          typename ARRAY_1D_W >
HPCREACT_HOST_DEVICE inline void
KineticReactions< REAL_TYPE,
                  INT_TYPE,
                  INDEX_TYPE,
                  ACTIVITY_MODEL,
                  LOGE_CONCENTRATION
                  >::computeReactionRatesElementary_impl( RealType const &, //temperature,
                                                          PARAMS_DATA const & params,
                                                          ARRAY_1D_TO_CONST const & activities,
                                                          RealType const waterActivity,
                                                          ARRAY_1D & reactionRates,
                                                          ARRAY_2D & dReactionRate_dActivities,
                                                          ARRAY_1D_W & dReactionRates_dWaterActivity )
{

  if constexpr( !CALCULATE_DERIVATIVES )
  {
    HPCREACT_UNUSED_VAR( dReactionRate_dActivities );
    HPCREACT_UNUSED_VAR( dReactionRates_dWaterActivity );
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
      RealType logProductActivityForward = 0.0;
      RealType logProductActivityReverse = 0.0;

      // build the products for the forward and reverse reaction rates
      for( IntType i = 0; i < PARAMS_DATA::numSpecies(); ++i )
      {

        RealType const s_ri = params.stoichiometricMatrix( r, i );

        if( s_ri < 0.0 )
        {
          logProductActivityForward += (-s_ri) * activities[i];
        }
        else if( s_ri > 0.0 )
        {
          logProductActivityReverse += s_ri * activities[i];
        }
      }

      // add water activity
      RealType const s_rw = params.waterStoichiometry( r );
      if( s_rw < 0.0 )
      {
        logProductActivityForward += (-s_rw) * waterActivity;
      }
      else if( s_rw > 0.0 )
      {
        logProductActivityReverse += s_rw * waterActivity;
      }

      reactionRates[r] = forwardRateConstant * exp( logProductActivityForward )
                         - reverseRateConstant * exp( logProductActivityReverse );

      if constexpr( CALCULATE_DERIVATIVES )
      {
        for( IntType i = 0; i < PARAMS_DATA::numSpecies(); ++i )
        {
          RealType const s_ri = params.stoichiometricMatrix( r, i );
          if( s_ri < 0.0 )
          {
            dReactionRate_dActivities[ r ][ i ] = forwardRateConstant * exp( logProductActivityForward ) * (-s_ri);
          }
          else if( s_ri > 0.0 )
          {
            dReactionRate_dActivities[ r ][ i ] = -reverseRateConstant * exp( logProductActivityReverse ) * s_ri;
          }
          else
          {
            dReactionRate_dActivities[ r ][ i ] = 0.0;
          }
        }

        if( s_rw < 0.0 )
        {
          dReactionRates_dWaterActivity[r] = forwardRateConstant * exp( logProductActivityForward ) * (-s_rw);
        }
        else if( s_rw > 0.0 )
        {
          dReactionRates_dWaterActivity[r] = -reverseRateConstant * exp( logProductActivityReverse ) * s_rw;
        }
        else
        {
          dReactionRates_dWaterActivity[r] = 0.0;
        }
      }
    }
    else
    {
      // variables used to build the product terms for the forward and reverse reaction rates
      RealType productActivityForward = 1.0;
      RealType productActivityReverse = 1.0;

      RealType dProductActivityForward_dActivities[PARAMS_DATA::numSpecies()];
      RealType dProductActivityReverse_dActivities[PARAMS_DATA::numSpecies()];
      for( IntType i = 0; i < PARAMS_DATA::numSpecies(); ++i )
      {
        dProductActivityForward_dActivities[i] = 1.0;
        dProductActivityReverse_dActivities[i] = 1.0;
      }

      // build the products for the forward and reverse reaction rates
      for( IntType i = 0; i < PARAMS_DATA::numSpecies(); ++i )
      {

        RealType const s_ri = params.stoichiometricMatrix( r, i );
        RealType const productTerm_i = activities[i] > 1e-100 ? pow( activities[i], fabs( s_ri ) ) : 0.0;

        if( s_ri < 0.0 )
        {
          productActivityForward *= productTerm_i;
        }
        else if( s_ri > 0.0 )
        {
          productActivityReverse *= productTerm_i;
        }

        if constexpr( CALCULATE_DERIVATIVES )
        {

          if( s_ri < 0.0 )
          {
            for( IntType j = 0; j < PARAMS_DATA::numSpecies(); ++j )
            {
              if( i==j )
              {
                dProductActivityForward_dActivities[j] *= -s_ri * pow( activities[i], -s_ri-1 );
                dProductActivityReverse_dActivities[j] = 0.0;
              }
              else
              {
                dProductActivityForward_dActivities[j] *= productTerm_i;
              }
            }
          }
          else if( s_ri > 0.0 )
          {
            for( IntType j = 0; j < PARAMS_DATA::numSpecies(); ++j )
            {
              if( i==j )
              {
                dProductActivityReverse_dActivities[j] *= s_ri * pow( activities[i], s_ri-1 );
                dProductActivityForward_dActivities[j] = 0.0;
              }
              else
              {
                dProductActivityReverse_dActivities[j] *= productTerm_i;
              }
            }
          }
          else
          {
            dProductActivityForward_dActivities[i] = 0.0;
            dProductActivityReverse_dActivities[i] = 0.0;
          }
        }
      }
      // add water activity
      RealType const s_rw = params.waterStoichiometry( r );
      if( s_rw < 0.0 )
      {
        RealType const productTerm_w = pow( waterActivity, -s_rw );
        productActivityForward *= productTerm_w;
        for( IntType j = 0; j < PARAMS_DATA::numSpecies(); ++j )
        {
          dProductActivityForward_dActivities[j] *= productTerm_w;
        }
      }
      else if( s_rw > 0.0 )
      {
        RealType const productTerm_w = pow( waterActivity, s_rw );
        productActivityReverse *= productTerm_w;
        for( IntType j = 0; j < PARAMS_DATA::numSpecies(); ++j )
        {
          dProductActivityReverse_dActivities[j] *= productTerm_w;
        }
      }

      reactionRates[r] = forwardRateConstant * productActivityForward - reverseRateConstant * productActivityReverse;

      if constexpr( CALCULATE_DERIVATIVES )
      {
        for( IntType i = 0; i < PARAMS_DATA::numSpecies(); ++i )
        {
          dReactionRate_dActivities[ r ][ i ] = forwardRateConstant * dProductActivityForward_dActivities[i] - reverseRateConstant * dProductActivityReverse_dActivities[i];
        }

        dReactionRates_dWaterActivity[r] = 0.0;
        if( waterActivity > 1e-100 )
        {
          if( s_rw < 0.0 )
          {
            dReactionRates_dWaterActivity[r] = forwardRateConstant * productActivityForward * (-s_rw) / waterActivity;
          }
          else if( s_rw > 0.0 )
          {
            dReactionRates_dWaterActivity[r] = -reverseRateConstant * productActivityReverse * s_rw / waterActivity;
          }
        }
      }
    } // end of if constexpr ( LOGE_CONCENTRATION )
  } // end of loop over reactions
}

template< typename REAL_TYPE,
          typename INT_TYPE,
          typename INDEX_TYPE,
          typename ACTIVITY_MODEL,
          bool LOGE_CONCENTRATION >
template< typename PARAMS_DATA,
          bool CALCULATE_DERIVATIVES,
          typename ARRAY_1D_TO_CONST,
          typename ARRAY_1D_SA,
          typename ARRAY_1D,
          typename ARRAY_2D,
          typename ARRAY_1D_W >
HPCREACT_HOST_DEVICE inline void
KineticReactions< REAL_TYPE,
                  INT_TYPE,
                  INDEX_TYPE,
                  ACTIVITY_MODEL,
                  LOGE_CONCENTRATION
                  >::computeReactionRatesAffinity_impl( RealType const &, //temperature,
                                                        PARAMS_DATA const & params,
                                                        ARRAY_1D_TO_CONST const & activities,
                                                        RealType const waterActivity,
                                                        ARRAY_1D_SA const & surfaceArea,
                                                        ARRAY_1D & reactionRates,
                                                        ARRAY_2D & dReactionRates_dActivities,
                                                        ARRAY_1D_W & dReactionRates_dWaterActivity )
{
  if constexpr( !CALCULATE_DERIVATIVES )
  {
    HPCREACT_UNUSED_VAR( dReactionRates_dActivities );
    HPCREACT_UNUSED_VAR( dReactionRates_dWaterActivity );
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
        dReactionRates_dActivities[ r ][ i ] = 0.0;
      }
    }

    // get/calculate the forward and reverse rate constants for this reaction
    RealType const rateConstant = params.rateConstantForward( r ); //* exp( -params.m_activationEnergy[r] / ( constants::R *
    // temperature ) );
    RealType const equilibriumConstant = params.equilibriumConstant( r );

    RealType quotient = 1.0;
    RealType const s_rw = params.waterStoichiometry( r );

    if constexpr( LOGE_CONCENTRATION )
    {
      RealType logQuotient = 0.0;
      // build the products for the forward and reverse reaction rates
      for( IntType i = 0; i < PARAMS_DATA::numSpecies(); ++i )
      {
        RealType const s_ri = params.stoichiometricMatrix( r, i );
        logQuotient += s_ri * activities[i];
      }
      // add water activity
      logQuotient += s_rw * waterActivity;
      quotient = exp( logQuotient );

      if constexpr( CALCULATE_DERIVATIVES )
      {
        for( IntType i = 0; i < PARAMS_DATA::numSpecies(); ++i )
        {
          RealType const s_ri = params.stoichiometricMatrix( r, i );
          dReactionRates_dActivities[ r ][ i ] = -rateConstant * surfaceArea[r] * s_ri * quotient / equilibriumConstant;
        }
        dReactionRates_dWaterActivity[r] = -rateConstant * surfaceArea[r] * s_rw * quotient / equilibriumConstant;
      } // end of if constexpr ( CALCULATE_DERIVATIVES )
    } // end of if constexpr ( LOGE_CONCENTRATION )
    else
    {
      for( IntType i = 0; i < PARAMS_DATA::numSpecies(); ++i )
      {

        RealType const s_ri = params.stoichiometricMatrix( r, i );
        if( s_ri > 0.0 || s_ri < 0.0 )
        {
          RealType const productTerm_i = activities[i] > 1e-100 ? pow( activities[i], s_ri ) : 0.0;
          quotient *= productTerm_i;
        }
      }
      // add water activity
      if( s_rw > 0.0 || s_rw < 0.0 )
      {
        quotient *= pow( waterActivity, s_rw );
      }

      if constexpr( CALCULATE_DERIVATIVES )
      {
        for( IntType i = 0; i < PARAMS_DATA::numSpecies(); ++i )
        {
          RealType const s_ri = params.stoichiometricMatrix( r, i );
          if( s_ri > 0.0 || s_ri < 0.0 )
          {
            dReactionRates_dActivities[ r ][ i ] = -rateConstant * surfaceArea[r] * s_ri * quotient / ( equilibriumConstant * activities[i] );
          }
          else
          {
            dReactionRates_dActivities[ r ][ i ] = 0.0;
          }
        }
        dReactionRates_dWaterActivity[r] = 0.0;
        if( waterActivity > 1e-100 )
        {
          dReactionRates_dWaterActivity[r] =
            -rateConstant * surfaceArea[r] * s_rw * quotient / ( equilibriumConstant * waterActivity );
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
          typename ACTIVITY_MODEL,
          bool LOGE_CONCENTRATION >
template< typename PARAMS_DATA,
          bool CALCULATE_DERIVATIVES,
          typename ARRAY_1D_TO_CONST,
          typename ARRAY_1D,
          typename ARRAY_2D,
          typename ARRAY_1D_W >
HPCREACT_HOST_DEVICE inline void
KineticReactions< REAL_TYPE,
                  INT_TYPE,
                  INDEX_TYPE,
                  ACTIVITY_MODEL,
                  LOGE_CONCENTRATION
                  >::computeSpeciesRates_impl( RealType const & temperature,
                                               PARAMS_DATA const & params,
                                               ARRAY_1D_TO_CONST const & activities,
                                               RealType const waterActivity,
                                               ARRAY_1D & speciesRates,
                                               ARRAY_2D & dSpeciesRates_dActivities,
                                               ARRAY_1D_W & dSpeciesRates_dWaterActivity )
{
  RealType reactionRates[PARAMS_DATA::numReactions()] = { 0.0 };
  RealType dReactionRates_dWaterActivity[PARAMS_DATA::numReactions()] = { 0.0 };
  CArrayWrapper< double, PARAMS_DATA::numReactions(), PARAMS_DATA::numSpecies() > dReactionRates_dActivities;

  if constexpr( !CALCULATE_DERIVATIVES )
  {
    HPCREACT_UNUSED_VAR( dSpeciesRates_dActivities );
    HPCREACT_UNUSED_VAR( dSpeciesRates_dWaterActivity );
  }

  computeReactionRatesElementary_impl< PARAMS_DATA, true >( temperature,
                                                            params,
                                                            activities,
                                                            waterActivity,
                                                            reactionRates,
                                                            dReactionRates_dActivities,
                                                            dReactionRates_dWaterActivity );

  for( IntType i = 0; i < PARAMS_DATA::numSpecies(); ++i )
  {
    speciesRates[i] = 0.0;
    if constexpr( CALCULATE_DERIVATIVES )
    {
      dSpeciesRates_dWaterActivity[i] = 0.0;
      for( IntType j = 0; j < PARAMS_DATA::numSpecies(); ++j )
      {
        dSpeciesRates_dActivities[ i ][ j ] = 0.0;
      }
    }
    for( IntType r=0; r<PARAMS_DATA::numReactions(); ++r )
    {
      RealType const s_ir = params.stoichiometricMatrix( r, i );
      speciesRates[i] += s_ir * reactionRates[r];
      if constexpr( CALCULATE_DERIVATIVES )
      {
        dSpeciesRates_dWaterActivity[i] += s_ir * dReactionRates_dWaterActivity[r];
        for( IntType j = 0; j < PARAMS_DATA::numSpecies(); ++j )
        {
          dSpeciesRates_dActivities[ i ][ j ] += s_ir * dReactionRates_dActivities[ r ][ j ];
        }
      }
    }
  }
}

template< typename REAL_TYPE,
          typename INT_TYPE,
          typename INDEX_TYPE,
          typename ACTIVITY_MODEL,
          bool LOGE_CONCENTRATION >
template< typename PARAMS_DATA,
          typename ARRAY_1D,
          typename ARRAY_1D_TO_CONST,
          typename ARRAY_2D >
HPCREACT_HOST_DEVICE inline void
KineticReactions< REAL_TYPE,
                  INT_TYPE,
                  INDEX_TYPE,
                  ACTIVITY_MODEL,
                  LOGE_CONCENTRATION >::timeStep( RealType const dt,
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
