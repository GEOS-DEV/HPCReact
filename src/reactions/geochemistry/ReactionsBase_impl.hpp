#include "ReactionsBase.hpp"

#include "common/macros.hpp"

namespace hpcReact
{

template< typename REAL_TYPE,
          typename REAL_DATA_ARRAY_1D_VIEW_TYPE,
          typename REAL_CONST_DATA_ARRAY_1D_VIEW_TYPE,
          typename INT_TYPE,
          typename INDEX_TYPE
          >
template< typename PARAMS_DATA >
HPCREACT_HOST_DEVICE inline
void ReactionsBase< REAL_TYPE,
                    REAL_DATA_ARRAY_1D_VIEW_TYPE,
                    REAL_CONST_DATA_ARRAY_1D_VIEW_TYPE,
                    INT_TYPE,
                    INDEX_TYPE >::computeLog10ActCoefBDotModel( REAL_TYPE const, //temperature,
                                                                REAL_TYPE const ionicStrength,
                                                                PARAMS_DATA const & params,
                                                                RealDataArrayView1d & log10PrimaryActCoeff,
                                                                RealDataArrayView1d & dLog10PrimaryActCoeff_dIonicStrength,
                                                                RealDataArrayView1d & log10SecActCoeff,
                                                                RealDataArrayView1d & dLog10SecActCoeff_dIonicStrength )
{
  // Compute log10(ActivityCoefficient) for basis and dependent species along with their
  // derivatives with respect to Ionic strength using the B-Dot Model
  // which is the same as the Extended Debye-Huckel model in GEOS.

  for( INDEX_TYPE i = 0; i < PARAMS_DATA::numPrimarySpecies; ++i )
  {
    log10PrimaryActCoeff[i] = params.m_WATEQBDot * ionicStrength - params.m_DebyeHuckelA * params.m_chargePrimary[i] * params.m_chargePrimary[i] * sqrt( ionicStrength ) /
                              (1.0 + params.m_ionSizePrimary[i] * params.m_DebyeHuckelB * sqrt( ionicStrength ));
    dLog10PrimaryActCoeff_dIonicStrength[i] = params.m_WATEQBDot - params.m_DebyeHuckelA * params.m_chargePrimary[i] * params.m_chargePrimary[i] *
                                              (0.5 / sqrt( ionicStrength ) / (1.0 + params.m_ionSizePrimary[i] * params.m_DebyeHuckelB * sqrt( ionicStrength )) - 0.5 * params.m_ionSizePrimary[i] *
                                               params.m_DebyeHuckelB /
                                               (1.0 + params.m_ionSizePrimary[i] * params.m_DebyeHuckelB * sqrt( ionicStrength )) /
                                               (1.0 + params.m_ionSizePrimary[i] * params.m_DebyeHuckelB * sqrt( ionicStrength )));
  }
  for( INDEX_TYPE i = 0; i < PARAMS_DATA::numSecondarySpecies; ++i )
  {
    log10SecActCoeff[i] = params.m_WATEQBDot * ionicStrength - params.m_DebyeHuckelA * params.m_chargeSec[i] * params.m_chargeSec[i] * sqrt( ionicStrength ) /
                          (1.0 + params.m_ionSizeSec[i] * params.m_DebyeHuckelB * sqrt( ionicStrength ));
    dLog10SecActCoeff_dIonicStrength[i] = params.m_WATEQBDot - params.m_DebyeHuckelA * params.m_chargeSec[i] * params.m_chargeSec[i] *
                                          (0.5 / sqrt( ionicStrength ) / (1.0 + params.m_ionSizeSec[i] * params.m_DebyeHuckelB * sqrt( ionicStrength )) - 0.5 * params.m_ionSizeSec[i] *
                                           params.m_DebyeHuckelB /
                                           (1.0 + params.m_ionSizeSec[i] * params.m_DebyeHuckelB * sqrt( ionicStrength )) /
                                           (1.0 + params.m_ionSizeSec[i] * params.m_DebyeHuckelB * sqrt( ionicStrength )));
  }
}


template< typename REAL_TYPE,
          typename REAL_DATA_ARRAY_1D_VIEW_TYPE,
          typename REAL_CONST_DATA_ARRAY_1D_VIEW_TYPE,
          typename INT_TYPE,
          typename INDEX_TYPE
          >
template< typename PARAMS_DATA >
HPCREACT_HOST_DEVICE inline
void ReactionsBase< REAL_TYPE,
                    REAL_DATA_ARRAY_1D_VIEW_TYPE,
                    REAL_CONST_DATA_ARRAY_1D_VIEW_TYPE,
                    INT_TYPE,
                    INDEX_TYPE >::computeIonicStrength( PARAMS_DATA const & params,
                                                        RealConstDataArrayView1d const & primarySpeciesConcentration,
                                                        RealConstDataArrayView1d const & secondarySpeciesConcentration,
                                                        REAL_TYPE & ionicStrength )
{
  //get ionic strength
  ionicStrength = 0.0;
  // Primary species
  for( INDEX_TYPE i = 0; i < PARAMS_DATA::numPrimarySpecies; ++i )
  {
    ionicStrength += 0.5 * params.m_chargePrimary[i] * params.m_chargePrimary[i] * primarySpeciesConcentration[i];
  }
  // Secondary species
  for( INDEX_TYPE j = 0; j < PARAMS_DATA::numSecondarySpecies; ++j )
  {
    ionicStrength += 0.5 * params.m_chargeSec[j] * params.m_chargeSec[j] * secondarySpeciesConcentration[j];
  }
}

}

#include "common/macrosCleanup.hpp"
