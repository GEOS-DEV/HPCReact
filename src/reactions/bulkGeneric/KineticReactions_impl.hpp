#include "KineticReactions.hpp"
#include "common/constants.hpp"
#include "common/macros.hpp"

#include <cmath>
#include <string>
#include <iostream>

#define OUTPUT_RATE_EQUATIONS 0
#if OUTPUT_RATE_EQUATIONS==1
  #include <sstream>
  #include <iomanip> 
  #define OUTPUT_RATE_FUNC(...) __VA_ARGS__
#else
  #define OUTPUT_RATE_FUNC(...)
#endif

namespace hpcReact
{
namespace bulkGeneric
{

template< typename REAL_TYPE,
          typename REAL_DATA_ARRAY_1D_VIEW_TYPE,
          typename REAL_CONST_DATA_ARRAY_1D_VIEW_TYPE,
          typename REAL_CONST_DATA_ARRAY_2D_VIEW_TYPE,
          typename INT_TYPE,
          typename INDEX_TYPE >
template< typename PARAMS_DATA, bool CALCULATE_DERIVATIVES >
HPCREACT_HOST_DEVICE inline void
KineticReactions< REAL_TYPE,
                  REAL_DATA_ARRAY_1D_VIEW_TYPE,
                  REAL_CONST_DATA_ARRAY_1D_VIEW_TYPE,
                  REAL_CONST_DATA_ARRAY_2D_VIEW_TYPE,
                  INT_TYPE,
                  INDEX_TYPE
                  >::computeReactionRates_impl( RealType const & ,//temperature,
                                                PARAMS_DATA const & params,
                                                RealConstDataArrayView1d const & primarySpeciesConcentration,
                                                RealDataArrayView1d & reactionRates,
                                                RealDataArrayView2d & reactionRatesDerivatives )
{

  OUTPUT_RATE_FUNC( 
    std::string const c[] = { "C02", "H2CO3", "Hp", "HCO3m", "CO32m", "Mg2SiO4", "Mg2p", "SiO2", "MgCO3", "H20" };
//    std::string const c[] = { "C02", "H2CO3", "H+", "HCO3-", "CO3^2-", "Mg2SiO4", "Mg2+", "SiO2", "MgCO3", "H20" };
    std::string reactionRatesString;
                  );
  OUTPUT_RATE_FUNC( std::string reactionRatesDerivativeString[PARAMS_DATA::numSpecies]; );

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

    RealType dProductConcForward_dC[PARAMS_DATA::numSpecies] = { 1.0 };
    RealType dProductConcReverse_dC[PARAMS_DATA::numSpecies] = { 1.0 };
    bool productConcForwardFlag = false;
    bool productConcReverseFlag = false;

    OUTPUT_RATE_FUNC( std::string productConcForwardString; std::string productConcReverseString; );
    OUTPUT_RATE_FUNC( std::string dProductConcForward_dCString[PARAMS_DATA::numSpecies]; std::string dProductConcReverse_dCString[PARAMS_DATA::numSpecies]; );

    for( IntType i = 0; i < PARAMS_DATA::numSpecies; ++i )
    {
      RealType const s_ri = params.stoichiometricMatrix( r, i );
      RealType const productTerm_i = pow( primarySpeciesConcentration[i], abs(s_ri) );
      
      if( s_ri < 0.0 )
      {
        productConcForward *= productTerm_i;
        productConcForwardFlag = true;
        OUTPUT_RATE_FUNC( { std::ostringstream oss;  oss << std::fixed << std::setprecision(0) << -s_ri; productConcForwardString += "*(" +c[i]+")" + "^" + oss.str() + " "; } );
      }
      else if( s_ri > 0.0 )
      {
        productConcReverse *= productTerm_i;
        productConcReverseFlag = true;
        OUTPUT_RATE_FUNC( { std::ostringstream oss;  oss << std::fixed << std::setprecision(0) << s_ri; productConcReverseString += "*(" +c[i]+")" + "^" + oss.str() + " "; } );
      }
      
      if constexpr ( CALCULATE_DERIVATIVES )
      {
        for( IntType j = 0; j < PARAMS_DATA::numSpecies; ++j )
        {
          if( i == j )
          {
            if( s_ri < 0.0 )
            {
              dProductConcForward_dC[j] *= -s_ri * productTerm_i / primarySpeciesConcentration[i];
              dProductConcReverse_dC[j] = 0.0;
              OUTPUT_RATE_FUNC( { std::ostringstream oss; oss <<std::fixed<<std::setprecision(0) << -s_ri-1; dProductConcForward_dCString[j] += "*"<<-s_ri<<"*(" +c[i]+")" + "^" + oss.str() + " "; } );
              OUTPUT_RATE_FUNC( dProductConcReverse_dCString[j] += "*0"; );
            }
            else if( s_ri > 0.0 )
            {
              dProductConcReverse_dC[j] *= s_ri * productTerm_i / primarySpeciesConcentration[i];
              dProductConcForward_dC[j] = 0.0;
              OUTPUT_RATE_FUNC( { std::ostringstream oss; oss <<std::fixed<<std::setprecision(0) << s_ri-1; dProductConcReverse_dCString[j] += "*"<<s_ri<<"*(" +c[i]+")" + "^" + oss.str() + " "; } );
              OUTPUT_RATE_FUNC( dProductConcForward_dCString[j] += "*0"; );
            }
            else
            {
              dProductConcForward_dC[j] =0.0;
              dProductConcReverse_dC[j] =0.0;
              OUTPUT_RATE_FUNC( { dProductConcForward_dCString[j] += "*0"; } );
              OUTPUT_RATE_FUNC( { dProductConcReverse_dCString[j] += "*0"; } );
            }
          }
          else
          {
            if( s_ri < 0.0 )
            {
              dProductConcForward_dC[j] *= productTerm_i;
              OUTPUT_RATE_FUNC( { std::ostringstream oss;  oss << std::fixed << std::setprecision(0) << -s_ri; dProductConcForward_dCString[j] += "*(" +c[i]+")" + "^" + oss.str() + " "; } );
            }
            else if( s_ri > 0.0 )
            {
              dProductConcReverse_dC[j] *= productTerm_i;
              OUTPUT_RATE_FUNC( { std::ostringstream oss; oss <<std::fixed<<std::setprecision(0) << s_ri; dProductConcReverse_dCString[j] += "*(" +c[i]+")" + "^" + oss.str() + " "; } );
            }
          }
        }
      }
    }

    
    OUTPUT_RATE_FUNC( reactionRatesString = "("; );
    if( productConcForwardFlag )
    {
      reactionRates[r] += forwardRateConstant * productConcForward;
      OUTPUT_RATE_FUNC( reactionRatesString += " k" + std::to_string(r+1) + "f" + productConcForwardString; );
    }
    if( productConcReverseFlag )
    {
      reactionRates[r] -= reverseRateConstant * productConcReverse;
      OUTPUT_RATE_FUNC( reactionRatesString += " - k" + std::to_string(r+1) + "r" + productConcReverseString; );
    }

    OUTPUT_RATE_FUNC( reactionRatesString += " )"; printf( "Reaction%d = %s\n", r, reactionRatesString.c_str() ); );

    if constexpr( CALCULATE_DERIVATIVES )
    {
      for( IntType i = 0; i < PARAMS_DATA::numSpecies; ++i )
      {
        OUTPUT_RATE_FUNC( reactionRatesDerivativeString[i] = "("; );
        if( productConcForwardFlag )
        {
          reactionRatesDerivatives( r, i ) += forwardRateConstant * dProductConcForward_dC[i];
          OUTPUT_RATE_FUNC( reactionRatesDerivativeString[i] = "k" + std::to_string(r+1) + "f" + dProductConcForward_dCString[i]; );
        }
        if( productConcReverseFlag )
        {
          reactionRatesDerivatives( r, i ) -= reverseRateConstant * dProductConcReverse_dC[i];
          OUTPUT_RATE_FUNC( reactionRatesDerivativeString[i] += " - k" + std::to_string(r+1) + "r" + dProductConcReverse_dCString[i]; );
        }

        OUTPUT_RATE_FUNC( printf( " dR%dd%s = %s\n", r, c[i].c_str(), reactionRatesDerivativeString[i].c_str() ); );

      }
    }

//    printf( " kf, kr, reactionRates[%2d] = % 8.4e, % 8.4e, % 8.4e\n", r, forwardRateConstant, reverseRateConstant, reactionRates[r] );
  }
}








// function to  the reaction rate. Includes impact of temperature, concentration, surface area, volume fraction and porosity
template< typename REAL_TYPE,
          typename REAL_DATA_ARRAY_1D_VIEW_TYPE,
          typename REAL_CONST_DATA_ARRAY_1D_VIEW_TYPE,
          typename REAL_DATA_ARRAY_2D_VIEW_TYPE,
          typename INT_TYPE,
          typename INDEX_TYPE >
template< typename PARAMS_DATA, bool CALCULATE_DERIVATIVES >
HPCREACT_HOST_DEVICE inline void
KineticReactions< REAL_TYPE,
                  REAL_DATA_ARRAY_1D_VIEW_TYPE,
                  REAL_CONST_DATA_ARRAY_1D_VIEW_TYPE,
                  REAL_DATA_ARRAY_2D_VIEW_TYPE,
                  INT_TYPE,
                  INDEX_TYPE
                  >::computeSpeciesRates_impl( RealType const & temperature,
                                           PARAMS_DATA const & params,
                                           RealConstDataArrayView1d const & primarySpeciesConcentration,
                                           RealDataArrayView1d & primarySpeciesRates,
                                           RealDataArrayView2d & primarySpeciesRatesDerivatives )
{
  RealType reactionRates[PARAMS_DATA::numReactions] = { 0.0 };
  RealDataArrayView2d reactionRatesDerivatives ;//= { { 0.0 } };

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
      if constexpr ( !CALCULATE_DERIVATIVES )
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
