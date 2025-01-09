#include "KineticReactions.hpp"
#include "common/macros.hpp"

namespace hpcReact
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
  /// 1. Create local vectors

  RealType log10PrimaryActCoeff[ PARAMS_DATA::numPrimarySpecies ];
  RealType log10SecActCoeff[ PARAMS_DATA::numSecondarySpecies ];

  RealType ionicStrength = 0.0;
  RealType dIonicStrength_dPrimaryConcentration[ PARAMS_DATA::numPrimarySpecies ];
  RealType dLog10PrimaryActCoeff_dIonicStrength[ PARAMS_DATA::numPrimarySpecies ];
  RealType dLog10SecActCoeff_dIonicStrength[ PARAMS_DATA::numSecondarySpecies ];

  /// 2. Compute activity coefficients
  computeIonicStrength( primarySpeciesConcentration,
                        secondarySpeciesConcentration,
                        ionicStrength );

  computeLog10ActCoefBDotModel( temperature,
                                ionicStrength,
                                log10PrimaryActCoeff,
                                dLog10PrimaryActCoeff_dIonicStrength,
                                log10SecActCoeff,
                                dLog10SecActCoeff_dIonicStrength );


  /// 3. Compute reaction rates
  for( int iRxn = 0; iRxn < PARAMS_DATA::numKineticReactions; iRxn++ )
  {
    RealType const saturationIndex = -params.m_log10EqConst[iRxn];

    for( int iPri = 0; iPri < PARAMS_DATA::numPrimarySpecies; ++iPri )
    {
      saturationIndex = saturationIndex 
                      + params.m_stoichMatrix[iRxn][iPri] * log10( primarySpeciesConcentration[iPri] ) 
                      + params.m_stoichMatrix[iRxn][iPri] * log10PrimaryActCoeff[iPri];
    }
    reactionRates[iRxn] = params.m_specificSurfaceArea * (1.0 - pow( 10, saturationIndex ) ) * params.m_reactionRateConstant[iRxn];
  }
}

} // namespace hpcReact

#include "common/macrosCleanup.hpp"
