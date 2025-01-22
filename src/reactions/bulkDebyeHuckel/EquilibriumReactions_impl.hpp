
#include "EquilibriumReactions.hpp"
#include "common/macros.hpp"

#include <iostream>
//#include "lapack.h"
namespace hpcReact
{
namespace bulkGeneric
{

template< typename REAL_TYPE,
          typename REAL_DATA_ARRAY_1D_VIEW_TYPE,
          typename REAL_CONST_DATA_ARRAY_1D_VIEW_TYPE,
          typename REAL_DATA_ARRAY_2D_VIEW_TYPE,
          typename REAL_CONST_DATA_ARRAY_2D_VIEW_TYPE,
          typename INT_TYPE,
          typename INDEX_TYPE >
template< typename PARAMS_DATA >
HPCREACT_HOST_DEVICE inline void
EquilibriumReactions< REAL_TYPE,
                      REAL_DATA_ARRAY_1D_VIEW_TYPE,
                      REAL_CONST_DATA_ARRAY_1D_VIEW_TYPE,
                      REAL_DATA_ARRAY_2D_VIEW_TYPE,
                      REAL_CONST_DATA_ARRAY_2D_VIEW_TYPE,
                      INT_TYPE,
                      INDEX_TYPE
                      >::updateConcentrations( RealType const temperature,
                                               PARAMS_DATA const &,
                                               RealConstDataArrayView1d & primarySpeciesTotalConcentration,
                                               RealDataArrayView1d & primarySpeciesConcentration,
                                               RealDataArrayView1d & secondarySpeciesConcentration )

{

  setInitialGuess( primarySpeciesTotalConcentration, primarySpeciesConcentration );

  REAL_TYPE matrix[ PARAMS_DATA::numPrimarySpecies ][ PARAMS_DATA::numPrimarySpecies ];
  REAL_TYPE rhs[ PARAMS_DATA::numPrimarySpecies ];
  REAL_TYPE solution[ PARAMS_DATA::numPrimarySpecies ];

  bool converged = false;
  for( int iteration = 0; iteration < m_maxNumIterations; iteration++ )
  {

    for( int i = 0; i< PARAMS_DATA::numPrimarySpecies; i++ )
    {
      rhs[i] = 0.0;
      solution[i] = 0.0;
      for( int j = 0; j< PARAMS_DATA::numPrimarySpecies; j++ )
      {
        matrix[ i ][ j ] = 0.0;
      }
    }

    assembleEquilibriumReactionSystem( temperature,
                                       primarySpeciesTotalConcentration,
                                       primarySpeciesConcentration,
                                       secondarySpeciesConcentration,
                                       matrix,
                                       rhs );

    REAL_TYPE const residualNorm = 0.0; //BlasLapackLA::vectorNorm2( rhs.toSliceConst() );

    if( residualNorm < m_newtonTol && iteration >= 1 )
    {
      converged = true;
      break;
    }

//    BlasLapackLA::solveLinearSystem( matrix, rhs, solution );

    updatePrimarySpeciesConcentrations( solution, primarySpeciesConcentration );
  }
  if( !converged )
  {
    std::cout<<"Equilibrium reactions did not converge."<<std::endl;
  }
  // GEOS_ERROR_IF( !converged, "Equilibrium reactions did not converge." );
}

template< typename REAL_TYPE,
          typename REAL_DATA_ARRAY_1D_VIEW_TYPE,
          typename REAL_CONST_DATA_ARRAY_1D_VIEW_TYPE,
          typename REAL_DATA_ARRAY_2D_VIEW_TYPE,
          typename REAL_CONST_DATA_ARRAY_2D_VIEW_TYPE,
          typename INT_TYPE,
          typename INDEX_TYPE >
template< typename PARAMS_DATA >
HPCREACT_HOST_DEVICE inline void
EquilibriumReactions< REAL_TYPE,
                      REAL_DATA_ARRAY_1D_VIEW_TYPE,
                      REAL_CONST_DATA_ARRAY_1D_VIEW_TYPE,
                      REAL_DATA_ARRAY_2D_VIEW_TYPE,
                      REAL_CONST_DATA_ARRAY_2D_VIEW_TYPE,
                      INT_TYPE,
                      INDEX_TYPE
                      >::assembleEquilibriumReactionSystem( RealType const temperature,
                                                            RealConstDataArrayView1d & primarySpeciesTotalConcentration,
                                                            RealConstDataArrayView1d & primarySpeciesConcentration,
                                                            RealDataArrayView1d & secondarySpeciesConcentration,
                                                            RealDataArrayView2d & matrix,
                                                            RealDataArrayView1d & rhs )
{
  constexpr IndexType numPrimarySpecies = PARAMS_DATA::numPrimarySpecies;
  constexpr IndexType numSecondarySpecies = PARAMS_DATA::numSecondarySpecies;
  RealType log10PrimaryActCoeff[ numPrimarySpecies ];

  RealType log10SecActCoeff[ numSecondarySpecies ];
  RealType dLog10SecConc_dLog10PrimaryConc[ numSecondarySpecies ][ numPrimarySpecies ];

  RealType totalConcentration[ numPrimarySpecies ];
  RealType dTotalConc_dLog10PrimaryConc[ numPrimarySpecies][ numPrimarySpecies ];

  RealType dIonicStrength_dPrimaryConcentration[ numPrimarySpecies ];
  RealType dLog10PrimaryActCoeff_dIonicStrength[ numPrimarySpecies ];
  RealType dLog10SecActCoeff_dIonicStrength[ numSecondarySpecies ];

  REAL_TYPE ionicStrength = 0.0;

  /// activity coefficients
  computeIonicStrength( primarySpeciesConcentration,
                        secondarySpeciesConcentration,
                        ionicStrength );

  computeLog10ActCoefBDotModel( temperature,
                                ionicStrength,
                                log10PrimaryActCoeff,
                                dLog10PrimaryActCoeff_dIonicStrength,
                                log10SecActCoeff,
                                dLog10SecActCoeff_dIonicStrength );

  computeSecondarySpeciesConcAndDerivative( temperature,
                                            log10PrimaryActCoeff,
                                            dLog10PrimaryActCoeff_dIonicStrength,
                                            log10SecActCoeff,
                                            dLog10SecActCoeff_dIonicStrength,
                                            primarySpeciesConcentration,
                                            secondarySpeciesConcentration,
                                            dLog10SecConc_dLog10PrimaryConc );

  computeTotalConcAndDerivative( temperature,
                                 primarySpeciesConcentration,
                                 secondarySpeciesConcentration,
                                 dLog10SecConc_dLog10PrimaryConc,
                                 totalConcentration,
                                 dTotalConc_dLog10PrimaryConc );

  for( int i=0; i<numPrimarySpecies; i++ )
  {
    rhs[i] = 1 - totalConcentration[i] / primarySpeciesTotalConcentration[i];
    rhs[i] = -rhs[i];
    for( int j=0; j<numPrimarySpecies; j++ )
    {
      matrix( i, j ) = -dTotalConc_dLog10PrimaryConc( i, j ) / primarySpeciesTotalConcentration[i];
    }
  }
}



// function to compute the derivative of the concentration of secondary species with respect to the concentration of the primary species.
template< typename REAL_TYPE,
          typename REAL_DATA_ARRAY_1D_VIEW_TYPE,
          typename REAL_CONST_DATA_ARRAY_1D_VIEW_TYPE,
          typename REAL_DATA_ARRAY_2D_VIEW_TYPE,
          typename REAL_CONST_DATA_ARRAY_2D_VIEW_TYPE,
          typename INT_TYPE,
          typename INDEX_TYPE >
template< typename PARAMS_DATA >
HPCREACT_HOST_DEVICE inline void
EquilibriumReactions< REAL_TYPE,
                      REAL_DATA_ARRAY_1D_VIEW_TYPE,
                      REAL_CONST_DATA_ARRAY_1D_VIEW_TYPE,
                      REAL_DATA_ARRAY_2D_VIEW_TYPE,
                      REAL_CONST_DATA_ARRAY_2D_VIEW_TYPE,
                      INT_TYPE,
                      INDEX_TYPE
                      >::computeSecondarySpeciesConcAndDerivative( RealType const temperature,
                                                                   PARAMS_DATA const & params,
                                                                   RealConstDataArrayView1d & log10PrimaryActCoeff,
                                                                   RealConstDataArrayView1d & dLog10PrimaryActCoeff_dIonicStrength,
                                                                   RealConstDataArrayView1d & log10SecActCoeff,
                                                                   RealConstDataArrayView1d & dLog10SecActCoeff_dIonicStrength,
                                                                   RealConstDataArrayView1d & primarySpeciesConcentration,
                                                                   RealDataArrayView1d & secondarySpeciesConcentration,
                                                                   RealDataArrayView2d & dLog10SecConc_dLog10PrimaryConc )
{
  GEOS_UNUSED_VAR( temperature );
  constexpr IndexType numPrimarySpecies = PARAMS_DATA::numPrimarySpecies;
  constexpr IndexType numSecondarySpecies = PARAMS_DATA::numSecondarySpecies;

  // Compute d(concentration of dependent species)/d(concentration of basis species)
  for( int iSec = 0; iSec < numSecondarySpecies; iSec++ )
  {
    REAL_TYPE log10SecConc = -params.m_log10EqConst[iSec] - log10SecActCoeff[iSec];

    for( int jPri = 0; jPri < numPrimarySpecies; jPri++ )
    {
      REAL_TYPE const dIonicStrength_dPrimaryConc = log( 10 ) * 0.5 * params.m_chargePrimary[jPri] * params.m_chargePrimary[jPri];

      log10SecConc += params.m_stoichMatrix[iSec][jPri] * ( log10( primarySpeciesConcentration[jPri] ) + log10PrimaryActCoeff[jPri] );
      dLog10SecConc_dLog10PrimaryConc[iSec][jPri] += params.m_stoichMatrix[iSec][jPri] - dLog10SecActCoeff_dIonicStrength[iSec] * primarySpeciesConcentration[jPri] *
                                                     dIonicStrength_dPrimaryConc;
      for( int kDerivative = 0; kDerivative < numPrimarySpecies; kDerivative++ )
      {
        // add contribution to the derivtive from all primary activity coefficients
        dLog10SecConc_dLog10PrimaryConc[iSec][jPri] += params.m_stoichMatrix[iSec][kDerivative] * dLog10PrimaryActCoeff_dIonicStrength[kDerivative] * primarySpeciesConcentration[jPri] *
                                                       dIonicStrength_dPrimaryConc;
      }

    }
    secondarySpeciesConcentration[iSec] = pow( 10, log10SecConc );
  }

}

template< typename REAL_TYPE,
          typename REAL_DATA_ARRAY_1D_VIEW_TYPE,
          typename REAL_CONST_DATA_ARRAY_1D_VIEW_TYPE,
          typename REAL_DATA_ARRAY_2D_VIEW_TYPE,
          typename REAL_CONST_DATA_ARRAY_2D_VIEW_TYPE,
          typename INT_TYPE,
          typename INDEX_TYPE >
template< typename PARAMS_DATA >
HPCREACT_HOST_DEVICE inline void
EquilibriumReactions< REAL_TYPE,
                      REAL_DATA_ARRAY_1D_VIEW_TYPE,
                      REAL_CONST_DATA_ARRAY_1D_VIEW_TYPE,
                      REAL_DATA_ARRAY_2D_VIEW_TYPE,
                      REAL_CONST_DATA_ARRAY_2D_VIEW_TYPE,
                      INT_TYPE,
                      INDEX_TYPE
                      >::computeTotalConcAndDerivative( RealType const temperature,
                                                        PARAMS_DATA const & params,
                                                        RealConstDataArrayView1d & primarySpeciesConcentration,
                                                        RealConstDataArrayView1d & secondarySpeciesConcentration,
                                                        RealConstDataArrayView2d & dLog10SecConc_dLog10PrimaryConc,
                                                        RealDataArrayView1d & totalConc,
                                                        RealDataArrayView2d & dTotalConc_dLog10PrimaryConc )
{
  GEOS_UNUSED_VAR( temperature );
  constexpr IndexType numPrimarySpecies = PARAMS_DATA::numPrimarySpecies;
  constexpr IndexType numSecondarySpecies = PARAMS_DATA::numSecondarySpecies;

  // This function computes the total concentration and its derivative with respect to log10(basis species concentrations).
  for( int iPri = 0; iPri < numPrimarySpecies; iPri++ )
  {
    totalConc[iPri] = primarySpeciesConcentration[iPri];
    // d(total concentration)/d(log10(concentration))
    dTotalConc_dLog10PrimaryConc[iPri][iPri] = log( 10.0 ) * primarySpeciesConcentration[iPri];
    // contribution from all dependent species
    for( int jSec = 0; jSec < numSecondarySpecies; jSec++ )
    {
      totalConc[iPri] += params.m_stoichMatrix[jSec][iPri] * secondarySpeciesConcentration[jSec];
      for( int kDerivative = 0; kDerivative < numPrimarySpecies; kDerivative++ )
      {
        // add contribution to the derivtive from dependent species via the chain rule
        dTotalConc_dLog10PrimaryConc[iPri][kDerivative] += params.m_stoichMatrix[jSec][iPri] * log( 10.0 ) *
                                                           secondarySpeciesConcentration[jSec] * dLog10SecConc_dLog10PrimaryConc[jSec][kDerivative];
      }
    }
  }
}

template< typename REAL_TYPE,
          typename REAL_DATA_ARRAY_1D_VIEW_TYPE,
          typename REAL_CONST_DATA_ARRAY_1D_VIEW_TYPE,
          typename REAL_DATA_ARRAY_2D_VIEW_TYPE,
          typename REAL_CONST_DATA_ARRAY_2D_VIEW_TYPE,
          typename INT_TYPE,
          typename INDEX_TYPE >
template< typename PARAMS_DATA >
HPCREACT_HOST_DEVICE inline void
EquilibriumReactions< REAL_TYPE,
                      REAL_DATA_ARRAY_1D_VIEW_TYPE,
                      REAL_CONST_DATA_ARRAY_1D_VIEW_TYPE,
                      REAL_DATA_ARRAY_2D_VIEW_TYPE,
                      REAL_CONST_DATA_ARRAY_2D_VIEW_TYPE,
                      INT_TYPE,
                      INDEX_TYPE
                      >::updatePrimarySpeciesConcentrations( RealConstDataArrayView1d solution,
                                                             RealDataArrayView1d & primarySpeciesConcentration )
{
  constexpr IndexType numPrimarySpecies = PARAMS_DATA::numPrimarySpecies;
  constexpr IndexType numSecondarySpecies = PARAMS_DATA::numSecondarySpecies;
  for( IndexType i = 0; i < numPrimarySpecies; i++ )
  {
    primarySpeciesConcentration[i] *= pow( 10, solution[i] );
  }
}


template< typename REAL_TYPE,
          typename REAL_DATA_ARRAY_1D_VIEW_TYPE,
          typename REAL_CONST_DATA_ARRAY_1D_VIEW_TYPE,
          typename REAL_DATA_ARRAY_2D_VIEW_TYPE,
          typename REAL_CONST_DATA_ARRAY_2D_VIEW_TYPE,
          typename INT_TYPE,
          typename INDEX_TYPE >
template< typename PARAMS_DATA >
HPCREACT_HOST_DEVICE inline void
EquilibriumReactions< REAL_TYPE,
                      REAL_DATA_ARRAY_1D_VIEW_TYPE,
                      REAL_CONST_DATA_ARRAY_1D_VIEW_TYPE,
                      REAL_DATA_ARRAY_2D_VIEW_TYPE,
                      REAL_CONST_DATA_ARRAY_2D_VIEW_TYPE,
                      INT_TYPE,
                      INDEX_TYPE
                      >::setInitialGuess( RealConstDataArrayView1d & primarySpeciesTotalConcentration,
                                          RealDataArrayView1d & primarySpeciesConcentration )
{
  for( IndexType i = 0; i < PARAMS_DATA::numPrimarySpecies; i++ )
  {
    primarySpeciesConcentration[i] = primarySpeciesTotalConcentration[i];
  }
  REAL_TYPE const hPlusConcentration = 2*primarySpeciesConcentration[2]-2*primarySpeciesConcentration[3]-primarySpeciesConcentration[4]+2*primarySpeciesConcentration[5]+primarySpeciesConcentration[6];
  if( hPlusConcentration < 0 )
  {
    primarySpeciesConcentration[0] = -hPlusConcentration;
  }
  else
  {
    primarySpeciesConcentration[0] = 1e-7;
  }

}

} // namespace bulkGeneric
} // namespace hpcReact

#include "common/macrosCleanup.hpp"
