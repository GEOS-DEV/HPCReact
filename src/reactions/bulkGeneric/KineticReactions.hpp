#pragma once

#include "common/macros.hpp"

namespace hpcReact
{
namespace bulkGeneric
{

template< typename REAL_TYPE,
          typename INT_TYPE,
          typename INDEX_TYPE >
class KineticReactions
{
public:

  using RealType = REAL_TYPE;
  // using RealDataArrayView1d = REAL_DATA_ARRAY_1D_VIEW_TYPE;
  // using RealConstDataArrayView1d = REAL_CONST_DATA_ARRAY_1D_VIEW_TYPE;
  // using RealDataArrayView2d = REAL_DATA_ARRAY_2D_VIEW_TYPE;
  using IntType = INT_TYPE;
  using IndexType = INDEX_TYPE;

  template< typename PARAMS_DATA,
            typename ARRAY_1D_TO_CONST,
            typename ARRAY_1D,
            typename ARRAY_2D >
  static HPCREACT_HOST_DEVICE inline void
  computeReactionRates( RealType const & temperature,
                        PARAMS_DATA const & params,
                        ARRAY_1D_TO_CONST const & primarySpeciesConcentration,
                        ARRAY_1D & reactionRates,
                        ARRAY_2D & reactionRatesDerivatives )
  { 
    computeReactionRates_impl< PARAMS_DATA, true >( temperature, 
                                                         params, 
                                                         primarySpeciesConcentration, 
                                                         reactionRates, 
                                                         reactionRatesDerivatives );
  }

  template< typename PARAMS_DATA,
            typename ARRAY_1D_TO_CONST,
            typename ARRAY_1D >
  static HPCREACT_HOST_DEVICE inline void
  computeReactionRates( RealType const & temperature,
                        PARAMS_DATA const & params,
                        ARRAY_1D_TO_CONST const & primarySpeciesConcentration,
                        ARRAY_1D & reactionRates)
  { 
    REAL_TYPE reactionRatesDerivatives[PARAMS_DATA::numReactions][PARAMS_DATA::numSpecies] = { 0.0 };
    computeReactionRates_impl< PARAMS_DATA, false >( temperature, 
                                                          params, 
                                                          primarySpeciesConcentration, 
                                                          reactionRates, 
                                                          reactionRatesDerivatives );
  }


  template< typename PARAMS_DATA,
            typename ARRAY_1D_TO_CONST,
            typename ARRAY_1D,
            typename ARRAY_2D >
  static HPCREACT_HOST_DEVICE inline void
  computeSpeciesRates( RealType const & temperature,
                       PARAMS_DATA const & params,
                       ARRAY_1D_TO_CONST const & primarySpeciesConcentration,
                       ARRAY_1D & speciesRates,
                       ARRAY_2D & primarySpeciesRatesDerivatives )
  {
    computeSpeciesRates_impl< PARAMS_DATA, true >( temperature, 
                                                        params, 
                                                        primarySpeciesConcentration, 
                                                        speciesRates, 
                                                        primarySpeciesRatesDerivatives );
  }

  template< typename PARAMS_DATA,
            typename ARRAY_1D_TO_CONST,
            typename ARRAY_1D >
  static HPCREACT_HOST_DEVICE inline void
  computeSpeciesRates( RealType const & temperature,
                       PARAMS_DATA const & params,
                       ARRAY_1D_TO_CONST const & primarySpeciesConcentration,
                       ARRAY_1D & speciesRates )
  {
    double primarySpeciesRatesDerivatives;
    computeSpeciesRates_impl< PARAMS_DATA, false >( temperature, 
                                                         params, 
                                                         primarySpeciesConcentration, 
                                                         speciesRates, 
                                                         primarySpeciesRatesDerivatives );
  }


private:
  template< typename PARAMS_DATA, 
            bool CALCULATE_DERIVATIVES, 
            typename ARRAY_1D_TO_CONST, 
            typename ARRAY_1D, 
            typename ARRAY_2D >
  static HPCREACT_HOST_DEVICE inline void
  computeReactionRates_impl( RealType const & temperature,
                        PARAMS_DATA const & params,
                        ARRAY_1D_TO_CONST const & primarySpeciesConcentration,
                        ARRAY_1D & reactionRates,
                        ARRAY_2D & reactionRatesDerivatives );

  template< typename PARAMS_DATA, 
            bool CALCULATE_DERIVATIVES, 
            typename ARRAY_1D_TO_CONST, 
            typename ARRAY_1D, 
            typename ARRAY_2D > 
  static HPCREACT_HOST_DEVICE inline void
  computeSpeciesRates_impl( RealType const & temperature,
                            PARAMS_DATA const & params,
                            ARRAY_1D_TO_CONST const & primarySpeciesConcentration,
                            ARRAY_1D & speciesRates,
                            ARRAY_2D & primarySpeciesRatesDerivatives );

};

} // namespace bulkGeneric
} // namespace hpcReact

#include "common/macrosCleanup.hpp"
