#pragma once

/**
 * @file EquilibriumReactions.hpp
 */

#include "ReactionsBase.hpp"


namespace hpcReact
{

template< typename REAL_TYPE, 
          template< int > typename REAL_DATA_ARRAY_VIEW_TYPE,
          template< int > typename REAL_CONST_DATA_ARRAY_VIEW_TYPE,
          typename INT_TYPE,
          typename INT_DATA_ARRAY_VIEW_TYPE,
          typename INT_CONST_DATA_ARRAY_VIEW_TYPE,
          typename INDEX_TYPE
        >
class EquilibriumReactions : public ReactionsBase< REAL_TYPE,
                                                   REAL_DATA_ARRAY_VIEW_TYPE,
                                                   REAL_CONST_DATA_ARRAY_VIEW_TYPE,
                                                   INT_TYPE,
                                                   INT_DATA_ARRAY_VIEW_TYPE,
                                                   INT_CONST_DATA_ARRAY_VIEW_TYPE,
                                                   INDEX_TYPE >
{
public:

  using Base = ReactionsBase< REAL_TYPE,
                              REAL_DATA_ARRAY_VIEW_TYPE,
                              REAL_CONST_DATA_ARRAY_VIEW_TYPE,
                              INT_TYPE,
                              INT_DATA_ARRAY_VIEW_TYPE,
                              INT_CONST_DATA_ARRAY_VIEW_TYPE,
                              INDEX_TYPE >;

  void updateConcentrations( real64 const temperature,
                             RealConstDataArrayView1d & primarySpeciesTotalConcentration,
                             RealDataArrayView1d & primarySpeciesContentration,
                             RealDataArrayView1d & secondarySpeciesConcentration ) const;
private:

  static constexpr integer m_maxNumIterations = MultiFluidConstants::maxNewtonIterations;
  static constexpr real64 m_newtonTol = 1e-6;

  void assembleEquilibriumReactionSystem( real64 const temperature,
                                          RealConstDataArrayView1d & primarySpeciesTotalConcentration,
                                          RealConstDataArrayView1d & primarySpeciesConcentration,
                                          RealDataArrayView1d & secondarySpeciesConcentration,
                                          RealDataArrayView2d & matrix,
                                          RealDataArrayView1d & rhs ) const;

  void computeSeondarySpeciesConcAndDerivative( real64 const temperature,
                                                RealConstDataArrayView1d & log10PrimaryActCoeff,
                                                RealConstDataArrayView1d & dLog10PrimaryActCoeff_dIonicStrength,
                                                RealConstDataArrayView1d & log10SecActCoeff,
                                                RealConstDataArrayView1d & dLog10SecActCoeff_dIonicStrength,
                                                RealConstDataArrayView1d & primarySpeciesConcentration,
                                                RealDataArrayView1d & secondarySpeciesConcentration,
                                                RealDataArrayView2d & dLog10SecConc_dLog10PrimaryConc ) const;

  void computeTotalConcAndDerivative( real64 const temperature,
                                      RealConstDataArrayView1d & primarySpeciesConcentration,
                                      RealConstDataArrayView1d & secondarySpeciesConcentration,
                                      RealConstDataArrayView2d & dLog10SecConc_dLog10PrimaryConc,
                                      RealDataArrayView1d & totalConc,
                                      RealDataArrayView2d & dTotalConc_dLog10PrimaryConc ) const;

  void updatePrimarySpeciesConcentrations( RealConstDataArrayView1d solution,
                                           RealDataArrayView1d & primarySpeciesConcentration ) const;

  void setInitialGuess( RealConstDataArrayView1d & primarySpeciesTotalConcentration,
                        RealDataArrayView1d & primarySpeciesConcentration ) const;


};


} // namespace hpcReact