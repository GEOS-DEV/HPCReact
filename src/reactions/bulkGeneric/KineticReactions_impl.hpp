#include "KineticReactions.hpp"
#include "common/constants.hpp"
#include "common/CArrayWrapper.hpp"
#include "common/macros.hpp"

#include <cmath>
#include <string>
#include <iostream>

namespace hpcReact
{
namespace bulkGeneric
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
                  >::computeReactionRates_impl( RealType const & ,//temperature,
                                                PARAMS_DATA const & params,
                                                ARRAY_1D_TO_CONST const & primarySpeciesConcentration,
                                                ARRAY_1D & reactionRates,
                                                ARRAY_2D & reactionRatesDerivatives )
{

  if constexpr ( !CALCULATE_DERIVATIVES )
  {
    HPCREACT_UNUSED_VAR( reactionRatesDerivatives );
  }


  for( IntType r=0; r<PARAMS_DATA::numReactions; ++r )
  {
    reactionRates[r] = 0.0;
    RealType const forwardRateConstant = params.rateConstantForward(r) ;//* exp( -params.m_activationEnergy[r] / ( constants::R * temperature ) );
    RealType const reverseRateConstant = params.rateConstantReverse(r) ;

    RealType productConcForward = 1.0;
    RealType productConcReverse = 1.0;

    RealType dProductConcForward_dC[PARAMS_DATA::numSpecies];
    RealType dProductConcReverse_dC[PARAMS_DATA::numSpecies];
    for( IntType i = 0; i < PARAMS_DATA::numSpecies; ++i )
    {
      dProductConcForward_dC[i] = 1.0;
      dProductConcReverse_dC[i] = 1.0;
    }

    for( IntType i = 0; i < PARAMS_DATA::numSpecies; ++i )
    {

      RealType const s_ri = params.stoichiometricMatrix( r, i );
      RealType const productTerm_i = primarySpeciesConcentration[i] > 1e-100 ? pow( primarySpeciesConcentration[i], abs(s_ri) ) : 0.0;
      
      if( s_ri < 0.0 )
      {
        productConcForward *= productTerm_i;
      }
      else if( s_ri > 0.0 )
      {
        productConcReverse *= productTerm_i;
      }
      
      if constexpr ( CALCULATE_DERIVATIVES )
      {

        if( s_ri < 0.0 )
        {
          for( IntType j = 0; j < PARAMS_DATA::numSpecies; ++j )
          {
            if( i==j )
            {
              dProductConcForward_dC[j] *= -s_ri * pow( primarySpeciesConcentration[i], -s_ri-1 );
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
          for( IntType j = 0; j < PARAMS_DATA::numSpecies; ++j )
          {  
            if( i==j )
            {
              dProductConcReverse_dC[j] *= s_ri * pow( primarySpeciesConcentration[i], s_ri-1 );
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
//    printf( " kf, kr, reactionRates[%2d] = % 8.4e, % 8.4e, % 8.4e\n", r, forwardRateConstant, reverseRateConstant, reactionRates[r] );


    if constexpr( CALCULATE_DERIVATIVES )
    {
      for( IntType i = 0; i < PARAMS_DATA::numSpecies; ++i )
      {
        reactionRatesDerivatives( r, i ) = forwardRateConstant * dProductConcForward_dC[i] - reverseRateConstant * dProductConcReverse_dC[i];
//        printf( " dR%ddC%d = % 3.1f * %3.2f - %3.1f * % 3.2f = % 3.2f \n", r, i, forwardRateConstant, dProductConcForward_dC[i], reverseRateConstant, dProductConcReverse_dC[i], reactionRatesDerivatives( r, i ) );
      }
    }


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
                                           ARRAY_1D_TO_CONST const & primarySpeciesConcentration,
                                           ARRAY_1D & primarySpeciesRates,
                                           ARRAY_2D & primarySpeciesRatesDerivatives )
{
  RealType reactionRates[PARAMS_DATA::numReactions] = { 0.0 };
  CArrayWrapper< double, PARAMS_DATA::numReactions, PARAMS_DATA::numSpecies > reactionRatesDerivatives;

  if constexpr ( !CALCULATE_DERIVATIVES )
  {
    HPCREACT_UNUSED_VAR(primarySpeciesRatesDerivatives);
  }

  computeReactionRates< PARAMS_DATA >( temperature, params, primarySpeciesConcentration, reactionRates, reactionRatesDerivatives );

  for( IntType i = 0; i < PARAMS_DATA::numSpecies; ++i )
  {
    primarySpeciesRates[i] = 0.0;
    if constexpr ( CALCULATE_DERIVATIVES )
    {
      for( IntType j = 0; j < PARAMS_DATA::numSpecies; ++j )
      {
        primarySpeciesRatesDerivatives( i, j ) = 0.0;
      }
    }
    for( IntType r=0; r<PARAMS_DATA::numReactions; ++r )
    {
      RealType const s_ir = params.stoichiometricMatrix( r, i );
      primarySpeciesRates[i] += s_ir * reactionRates[r];
      if constexpr ( CALCULATE_DERIVATIVES )
      {
        for( IntType j = 0; j < PARAMS_DATA::numSpecies; ++j )
        {
          primarySpeciesRatesDerivatives( i, j ) += s_ir * reactionRatesDerivatives( r, j );
        }
      }
    }
  }
}

} // namespace bulkGeneric
} // namespace hpcReact

#include "common/macrosCleanup.hpp"
