#pragma once

#include "common/macros.hpp"

namespace hpcReact
{
namespace bulkGeneric
{

template< typename REAL_TYPE,
          typename REAL_DATA_ARRAY_1D_VIEW_TYPE,
          typename REAL_CONST_DATA_ARRAY_1D_VIEW_TYPE,
          typename REAL_DATA_ARRAY_2D_VIEW_TYPE,
          typename INT_TYPE,
          typename INDEX_TYPE >
class KineticReactions
{
public:

  using RealType = REAL_TYPE;
  using RealDataArrayView1d = REAL_DATA_ARRAY_1D_VIEW_TYPE;
  using RealConstDataArrayView1d = REAL_CONST_DATA_ARRAY_1D_VIEW_TYPE;
  using RealDataArrayView2d = REAL_DATA_ARRAY_2D_VIEW_TYPE;
  using IntType = INT_TYPE;
  using IndexType = INDEX_TYPE;

  template< typename PARAMS_DATA >
  static HPCREACT_HOST_DEVICE inline void
  computeReactionRates( RealType const & temperature,
                        PARAMS_DATA const & params,
                        RealConstDataArrayView1d const & primarySpeciesConcentration,
                        RealDataArrayView1d & reactionRates,
                        RealDataArrayView2d & reactionRatesDerivatives )
  { 
    computeReactionRates_impl< PARAMS_DATA, true >( temperature, 
                                                         params, 
                                                         primarySpeciesConcentration, 
                                                         reactionRates, 
                                                         reactionRatesDerivatives );
  }

  template< typename PARAMS_DATA >
  static HPCREACT_HOST_DEVICE inline void
  computeReactionRates( RealType const & temperature,
                        PARAMS_DATA const & params,
                        RealConstDataArrayView1d const & primarySpeciesConcentration,
                        RealDataArrayView1d & reactionRates)
  { 
    RealDataArrayView2d reactionRatesDerivatives;
    computeReactionRates_impl< PARAMS_DATA, false >( temperature, 
                                                          params, 
                                                          primarySpeciesConcentration, 
                                                          reactionRates, 
                                                          reactionRatesDerivatives );
  }


  template< typename PARAMS_DATA >
  static HPCREACT_HOST_DEVICE inline void
  computeSpeciesRates( RealType const & temperature,
                       PARAMS_DATA const & params,
                       RealConstDataArrayView1d const & primarySpeciesConcentration,
                       RealDataArrayView1d & speciesRates,
                       RealDataArrayView2d & primarySpeciesRatesDerivatives )
  {
    computeSpeciesRates_impl< PARAMS_DATA, true >( temperature, 
                                                        params, 
                                                        primarySpeciesConcentration, 
                                                        speciesRates, 
                                                        primarySpeciesRatesDerivatives );
  }

  template< typename PARAMS_DATA >
  static HPCREACT_HOST_DEVICE inline void
  computeSpeciesRates( RealType const & temperature,
                       PARAMS_DATA const & params,
                       RealConstDataArrayView1d const & primarySpeciesConcentration,
                       RealDataArrayView1d & speciesRates )
  {
    RealDataArrayView2d primarySpeciesRatesDerivatives;
    computeSpeciesRates_impl< PARAMS_DATA, false >( temperature, 
                                                         params, 
                                                         primarySpeciesConcentration, 
                                                         speciesRates, 
                                                         primarySpeciesRatesDerivatives );
  }


private:
  template< typename PARAMS_DATA, bool CALCULATE_DERIVATIVES >
  static HPCREACT_HOST_DEVICE inline void
  computeReactionRates_impl( RealType const & temperature,
                        PARAMS_DATA const & params,
                        RealConstDataArrayView1d const & primarySpeciesConcentration,
                        RealDataArrayView1d & reactionRates,
                        RealDataArrayView2d & reactionRatesDerivatives );

  template< typename PARAMS_DATA, bool CALCULATE_DERIVATIVES >
  static HPCREACT_HOST_DEVICE inline void
  computeSpeciesRates_impl( RealType const & temperature,
                       PARAMS_DATA const & params,
                       RealConstDataArrayView1d const & primarySpeciesConcentration,
                       RealDataArrayView1d & speciesRates,
                       RealDataArrayView2d & primarySpeciesRatesDerivatives );

};

} // namespace bulkGeneric
} // namespace hpcReact

#include "common/macrosCleanup.hpp"
