#pragma once

#include "common/macros.hpp"
#include "common/CArrayWrapper.hpp"
#include "common/DirectSystemSolve.hpp"
#include "common/printers.hpp"

#include <iostream>


namespace hpcReact
{
namespace bulkGeneric
{

template< typename REAL_TYPE,
          typename INT_TYPE,
          typename INDEX_TYPE >
class EquilibriumReactions
{
public:
  using RealType = REAL_TYPE;
  using IntType = INT_TYPE;
  using IndexType = INDEX_TYPE;
//  constexpr static int formulation = FORMULATION;

  template< typename PARAMS_DATA,
            typename ARRAY_1D,
            typename ARRAY_1D_TO_CONST,
            typename ARRAY_1D_TO_CONST2,
            typename ARRAY_2D >
  static HPCREACT_HOST_DEVICE void
  computeResidualAndJacobian( RealType const & temperature,
                              PARAMS_DATA const & params,
                              ARRAY_1D_TO_CONST const & speciesConcentration0,
                              ARRAY_1D_TO_CONST2 const & xi,
                              ARRAY_1D & residual,
                              ARRAY_2D & jacobian );


  template< typename PARAMS_DATA,
            typename ARRAY_1D,
            typename ARRAY_1D_TO_CONST >
  static HPCREACT_HOST_DEVICE
  void
  enforceEquilibrium_Extents( RealType const & temperature,
                      PARAMS_DATA const & params,
                      ARRAY_1D_TO_CONST const & speciesConcentration0,
                      ARRAY_1D & speciesConcentration );

  template< typename PARAMS_DATA,
            typename ARRAY_1D,
            typename ARRAY_1D_TO_CONST >
  static HPCREACT_HOST_DEVICE
  void
  enforceEquilibrium_Aggregate( RealType const & temperature,
                      PARAMS_DATA const & params,
                      ARRAY_1D_TO_CONST const & speciesConcentration0,
                      ARRAY_1D & speciesConcentration );

private:
  template< typename PARAMS_DATA,
            typename ARRAY_1D,
            typename ARRAY_1D_TO_CONST,
            typename ARRAY_1D_TO_CONST2,
            typename ARRAY_2D >
  static HPCREACT_HOST_DEVICE void
  computeResidualAndJacobianReactionExtents( RealType const & temperature,
                                             PARAMS_DATA const & params,
                                             ARRAY_1D_TO_CONST const & speciesConcentration0,
                                             ARRAY_1D_TO_CONST2 const & xi,
                                             ARRAY_1D & residual,
                                             ARRAY_2D & jacobian );

  template< typename PARAMS_DATA,
            typename ARRAY_1D,
            typename ARRAY_1D_TO_CONST,
            typename ARRAY_2D >
  static HPCREACT_HOST_DEVICE void
  computeResidualAndJacobianAggregatePrimaryConcentrations( RealType const & temperature,
                                                            PARAMS_DATA const & params,
                                                            ARRAY_1D_TO_CONST const & targetAggregatePrimaryConcentrations,
                                                            ARRAY_1D_TO_CONST const & logPrimarySpeciesConcentration,
                                                            ARRAY_1D & residual,
                                                            ARRAY_2D & jacobian );
};




} // namespace bulkGeneric
} // namespace hpcReact

#if !defined(__INTELLISENSE__)
#include "EquilibriumReactionsAggregatePrimaryConcentration_impl.hpp"
#include "EquilibriumReactionsReactionExtents_impl.hpp"
#endif
