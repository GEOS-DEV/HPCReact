#pragma once

/**
 * @file EquilibriumReactions.hpp
 */

#include "ReactionsBase.hpp"
#include "common/macros.hpp"

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
class EquilibriumReactions : public ReactionsBase< REAL_TYPE,
                                                   REAL_DATA_ARRAY_1D_VIEW_TYPE,
                                                   REAL_CONST_DATA_ARRAY_1D_VIEW_TYPE,
                                                   INT_TYPE,
                                                   INDEX_TYPE >
{
public:

  using Base = ReactionsBase< REAL_TYPE,
                              REAL_DATA_ARRAY_1D_VIEW_TYPE,
                              REAL_CONST_DATA_ARRAY_1D_VIEW_TYPE,
                              INT_TYPE,
                              INDEX_TYPE >;

  using typename Base::RealType;
  using typename Base::RealDataArrayView1d;
  using typename Base::RealConstDataArrayView1d;
  using typename Base::IntType;
  using typename Base::IndexType;

  using RealDataArrayView2d = REAL_DATA_ARRAY_2D_VIEW_TYPE;
  using RealConstDataArrayView2d = REAL_CONST_DATA_ARRAY_2D_VIEW_TYPE;

  template< typename PARAMS_DATA >
  static void updateConcentrations( RealType const temperature,
                                    PARAMS_DATA const & params,
                                    RealConstDataArrayView1d & primarySpeciesTotalConcentration,
                                    RealDataArrayView1d & primarySpeciesContentration,
                                    RealDataArrayView1d & secondarySpeciesConcentration );
private:

  static constexpr IntType m_maxNumIterations = 10;
  static constexpr RealType m_newtonTol = 1e-6;

  template< typename PARAMS_DATA >
  static HPCREACT_HOST_DEVICE inline void
  setInitialGuess( RealConstDataArrayView1d & primarySpeciesTotalConcentration,
                   RealDataArrayView1d & primarySpeciesConcentration );

  template< typename PARAMS_DATA >
  static HPCREACT_HOST_DEVICE inline void
  assembleEquilibriumReactionSystem( RealType const temperature,
                                     RealConstDataArrayView1d & primarySpeciesTotalConcentration,
                                     RealConstDataArrayView1d & primarySpeciesConcentration,
                                     RealDataArrayView1d & secondarySpeciesConcentration,
                                     RealDataArrayView2d & matrix,
                                     RealDataArrayView1d & rhs );

  template< typename PARAMS_DATA >
  static HPCREACT_HOST_DEVICE inline void
  computeSecondarySpeciesConcAndDerivative( RealType const temperature,
                                            PARAMS_DATA const & params,
                                            RealConstDataArrayView1d & log10PrimaryActCoeff,
                                            RealConstDataArrayView1d & dLog10PrimaryActCoeff_dIonicStrength,
                                            RealConstDataArrayView1d & log10SecActCoeff,
                                            RealConstDataArrayView1d & dLog10SecActCoeff_dIonicStrength,
                                            RealConstDataArrayView1d & primarySpeciesConcentration,
                                            RealDataArrayView1d & secondarySpeciesConcentration,
                                            RealDataArrayView2d & dLog10SecConc_dLog10PrimaryConc );

  template< typename PARAMS_DATA >
  static HPCREACT_HOST_DEVICE inline void
  computeTotalConcAndDerivative( RealType const temperature,
                                 PARAMS_DATA const & params,
                                 RealConstDataArrayView1d & primarySpeciesConcentration,
                                 RealConstDataArrayView1d & secondarySpeciesConcentration,
                                 RealConstDataArrayView2d & dLog10SecConc_dLog10PrimaryConc,
                                 RealDataArrayView1d & totalConc,
                                 RealDataArrayView2d & dTotalConc_dLog10PrimaryConc );

  template< typename PARAMS_DATA >
  static HPCREACT_HOST_DEVICE inline void
  updatePrimarySpeciesConcentrations( RealConstDataArrayView1d solution,
                                      RealDataArrayView1d & primarySpeciesConcentration );



};

} // namespace bulkGeneric
} // namespace hpcReact

#include "common/macrosCleanup.hpp"
