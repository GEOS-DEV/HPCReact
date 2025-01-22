#include "KineticReactions.hpp"
#include "common/macros.hpp"

namespace hpcReact
{
namespace bulkDebyeHuckel
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
                                           RealConstDataArrayView1d & primarySpeciesConcentration,
                                           RealConstDataArrayView1d & secondarySpeciesConcentration,
                                           RealDataArrayView1d & reactionRates )
{
  for( IntType i = 0; i < PARAMS_DATA::numPrimarySpecies; ++i )
  {
    reactionRates[i] = 0.0;
    for( IntType r=0; r<PARAMS_DATA::numKineticReactions; ++r )
    {
      RealType const forwardRateConstant = params.m_reactionRateConstant[r] * exp( -params.m_activationEnergy[r] / ( constants::R * temperature ) );
      RealType const reverseRateConstant = params.equilibriumConstant[r] / forwardRateConstant;
      RealType const s_ir = params.m_stoichiometricMatrix[r][i];

      RealType productConcPlus = 1.0;
      RealType productConcMinus = 1.0;
      bool productConcPlusFlag = false;
      bool productConcMinusFlag = false;
      for( IntType j = 0; j < PARAMS_DATA::numPrimarySpecies; ++j )
      {
        RealType const s_jr = params.m_stoichiometricMatrix[r][j];
        if( s_jr < 0.0 )
        {
          productConcPlus *= pow( primarySpeciesConcentration[j], -s_jr );
          productConcPlusFlag = true;
        }
        else if( s_jr > 0.0 )
        {
          productConcMinus *= pow( primarySpeciesConcentration[j], s_jr );
          productConcMinusFlag = true;
        }
      }

      RealType reactionRate_n = 0.0;
      if( productConcPlusFlag )
      {
        reactionRate_n += forwardRateConstant * productConcPlus;
      }
      if( productConcMinusFlag )
      {
        reactionRate_n -= reverseRateConstant * productConcMinus;
      }
      reactionRates[i] += s_ir * reactionRate_n;
    }
  }
}

} // namespace bulkGeneric
} // namespace hpcReact

#include "common/macrosCleanup.hpp"
