#include "KineticReactions.hpp"
#include "common/constants.hpp"
#include "common/macros.hpp"

#include <cmath>
#include <string>
#include <iostream>

#define OUTPUT_RATE_EQUATIONS 0

namespace hpcReact
{
namespace bulkGeneric
{

// function to  the reaction rate. Includes impact of temperature, concentration, surface area, volume fraction and porosity
template< typename REAL_TYPE,
          typename REAL_DATA_ARRAY_1D_VIEW_TYPE,
          typename REAL_CONST_DATA_ARRAY_1D_VIEW_TYPE,
          typename INT_TYPE,
          typename INDEX_TYPE >
template< typename PARAMS_DATA >
HPCREACT_HOST_DEVICE inline void
KineticReactions< REAL_TYPE,
                  REAL_DATA_ARRAY_1D_VIEW_TYPE,
                  REAL_CONST_DATA_ARRAY_1D_VIEW_TYPE,
                  INT_TYPE,
                  INDEX_TYPE
                  >::computeReactionRates( RealType const & temperature,
                                           PARAMS_DATA const & params,
                                           RealConstDataArrayView1d const & primarySpeciesConcentration,
                                           RealDataArrayView1d & reactionRates )
{
#if OUTPUT_RATE_EQUATIONS==1
  std::string const c[] = { "A", "B", "C", "D", "E" };
#endif

  for( IntType i = 0; i < PARAMS_DATA::numPrimarySpecies; ++i )
  {
#if OUTPUT_RATE_EQUATIONS==1
    std::string rateString;
#endif
    reactionRates[i] = 0.0;
    for( IntType r=0; r<PARAMS_DATA::numKineticReactions; ++r )
    {
      RealType const forwardRateConstant = params.m_rateConstant[r] * exp( -params.m_activationEnergy[r] / ( constants::R * temperature ) );
      RealType const reverseRateConstant = forwardRateConstant / params.m_equilibriumConstant[r] ;

//      printf( "forwardRateConstant = %8.4e, reverseRateConstant = %8.4e\n", forwardRateConstant, reverseRateConstant );
      RealType const s_ir = params.m_stoichiometricMatrix[r][i];

      RealType productConcPlus = 1.0;
      RealType productConcMinus = 1.0;
      bool productConcPlusFlag = false;
      bool productConcMinusFlag = false;
#if OUTPUT_RATE_EQUATIONS==1
      std::string productConcPlusString;
      std::string productConcMinusString;
#endif
      for( IntType j = 0; j < PARAMS_DATA::numPrimarySpecies; ++j )
      {
        RealType const s_jr = params.m_stoichiometricMatrix[r][j];
        if( s_jr < 0.0 )
        {
          productConcPlus *= pow( primarySpeciesConcentration[j], -s_jr );
#if OUTPUT_RATE_EQUATIONS==1
          productConcPlusString += " c(" +c[j]+")" ;//+ "^" + std::to_string(-s_jr) + " ";
#endif
          productConcPlusFlag = true;
        }
        else if( s_jr > 0.0 )
        {
          productConcMinus *= pow( primarySpeciesConcentration[j], s_jr );
#if OUTPUT_RATE_EQUATIONS==1
          productConcMinusString += " c("+c[j]+")" ;//+ "^" + std::to_string(s_jr) + " ";
#endif
          productConcMinusFlag = true;
        }
      }

      RealType reactionRate_n = 0.0;
#if OUTPUT_RATE_EQUATIONS==1
      std::string reactionRateString = "{ ";
#endif
      if( productConcPlusFlag )
      {
        reactionRate_n += forwardRateConstant * productConcPlus;
#if OUTPUT_RATE_EQUATIONS==1
        reactionRateString += " k_" + std::to_string(r+1) + "f" + productConcPlusString;
#endif
      }
      if( productConcMinusFlag )
      {
        reactionRate_n -= reverseRateConstant * productConcMinus;
#if OUTPUT_RATE_EQUATIONS==1
        reactionRateString += " - k_" + std::to_string(r+1) + "r" + productConcMinusString;
#endif
      }

#if OUTPUT_RATE_EQUATIONS==1
      reactionRateString += " }";
      if( s_ir > 0.1 || s_ir < -0.1 )
      {
        std::string s_irString;
        s_irString.resize(10);
        snprintf(s_irString.data(), 5, "%+.1f", s_ir);
        rateString += s_irString + " " + reactionRateString;
      }
#endif
      reactionRates[i] += s_ir * reactionRate_n;

    }
#if OUTPUT_RATE_EQUATIONS==1
    std::cout<<rateString<<std::endl;
#endif
  }
}

} // namespace bulkGeneric
} // namespace hpcReact

#include "common/macrosCleanup.hpp"
