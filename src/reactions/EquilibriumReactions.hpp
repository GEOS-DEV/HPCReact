#pragma once

/**
 * @file EquilibriumReactions.hpp
 */

#include "ReactionsBase.hpp"


namespace hpcReact
{

template< typename REAL_DATA_ARRAY_VIEW_TYPE,
          typename INTEGER_DATA_ARRAY_VIEW_TYPE,
          typename REAL_TYPE, 
          typename INT_TYPE,
          typename INDEX_TYPE
        >
class EquilibriumReactions : public ReactionsBase< REAL_DATA_ARRAY_VIEW_TYPE,
                                                   INTEGER_DATA_ARRAY_VIEW_TYPE,
                                                   REAL_TYPE, 
                                                   INT_TYPE,
                                                   INDEX_TYPE >
{
public:

  using Base = ReactionsBase< REAL_DATA_ARRAY_VIEW_TYPE,
                              INTEGER_DATA_ARRAY_VIEW_TYPE,
                              REAL_TYPE, 
                              INT_TYPE,
                              INDEX_TYPE >;

  using ParamsData = typename Base::ParamsData;


  void updateConcentrations( real64 const temperature,
                             arraySlice1d< real64 const, compflow::USD_COMP - 1 > const & primarySpeciesTotalConcentration,
                             arraySlice1d< real64, compflow::USD_COMP - 1 > const & primarySpeciesContentration,
                             arraySlice1d< real64, compflow::USD_COMP - 1 > const & secondarySpeciesConcentration ) const;
private:

  static constexpr integer m_maxNumIterations = MultiFluidConstants::maxNewtonIterations;
  static constexpr real64 m_newtonTol = 1e-6;

  void assembleEquilibriumReactionSystem( real64 const temperature,
                                          arraySlice1d< real64 const, compflow::USD_COMP - 1 > const & primarySpeciesTotalConcentration,
                                          arraySlice1d< real64 const, compflow::USD_COMP - 1 > const & primarySpeciesConcentration,
                                          arraySlice1d< real64, compflow::USD_COMP - 1 > const & secondarySpeciesConcentration,
                                          arraySlice2d< real64 > const & matrix,
                                          arraySlice1d< real64 > const & rhs ) const;

  void computeSeondarySpeciesConcAndDerivative( real64 const temperature,
                                                arraySlice1d< real64 const > const & log10PrimaryActCoeff,
                                                arraySlice1d< real64 const > const & dLog10PrimaryActCoeff_dIonicStrength,
                                                arraySlice1d< real64 const > const & log10SecActCoeff,
                                                arraySlice1d< real64 const > const & dLog10SecActCoeff_dIonicStrength,
                                                arraySlice1d< real64 const, compflow::USD_COMP - 1 > const & primarySpeciesConcentration,
                                                arraySlice1d< real64, compflow::USD_COMP - 1 > const & secondarySpeciesConcentration,
                                                arraySlice2d< real64 > const & dLog10SecConc_dLog10PrimaryConc ) const;

  void computeTotalConcAndDerivative( real64 const temperature,
                                      arraySlice1d< real64 const, compflow::USD_COMP - 1 > const & primarySpeciesConcentration,
                                      arraySlice1d< real64 const, compflow::USD_COMP - 1 > const & secondarySpeciesConcentration,
                                      arraySlice2d< real64 const > const & dLog10SecConc_dLog10PrimaryConc,
                                      arraySlice1d< real64 > const & totalConc,
                                      arraySlice2d< real64 > const & dTotalConc_dLog10PrimaryConc ) const;

  void updatePrimarySpeciesConcentrations( arraySlice1d< real64 const > const solution,
                                           arraySlice1d< real64, compflow::USD_COMP - 1 > const & primarySpeciesConcentration ) const;

  void setInitialGuess( arraySlice1d< real64 const, compflow::USD_COMP - 1 > const & primarySpeciesTotalConcentration,
                        arraySlice1d< real64, compflow::USD_COMP - 1 > const & primarySpeciesConcentration ) const;


};


} // namespace hpcReact