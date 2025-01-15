#include "KineticReactions.hpp"
#include "common/macros.hpp"

namespace hpcReact
{
namespace liquidElectrolyteBatteryChemistry
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
  for( int iRxn = 0; iRxn < PARAMS_DATA::numKineticReactions; iRxn++ )
  {
    RealType const forwardRateConstant = params.m_reactionRateConstant[iRxn] * exp( -params.m_activationEnergy[iRxn] / ( constants::R * temperature ) );
    RealType const reverseRateConstant = params.equilibriumConstant[iRxn] / forwardRateConstant;

    for( int iPri = 0; iPri < PARAMS_DATA::numPrimarySpecies; ++iPri )
    {
      RealType const stoichiometricCoefficient = params.m_stoichiometricMatrix[iRxn][iPri];
      RealType const primarySpeciesConcentration = primarySpeciesConcentration[iPri];
      RealType const secondarySpeciesConcentration = secondarySpeciesConcentration[iPri];

      if( stoichiometricCoefficient < 0.0 )
      {
        reactionRates[iPri] = reactionRates[iRxn] + stoichiometricCoefficient * forwardRateConstant * primarySpeciesConcentration;
      }
      else if( stoichiometricCoefficient > 0.0 )
      {
        reactionRates[iPri] = reactionRates[iRxn] - stoichiometricCoefficient * reverseRateConstant * primarySpeciesConcentration;
      }
    }
  }
}

} // namespace geochemistry
} // namespace hpcReact

#include "common/macrosCleanup.hpp"
