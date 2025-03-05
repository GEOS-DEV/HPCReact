#pragma once

#include "common/macros.hpp"
#include "common/CArrayWrapper.hpp"
#include "common/DirectSystemSolve.hpp"


namespace hpcReact
{
namespace bulkGeneric
{

template< typename REAL_TYPE,
          typename INT_TYPE,
          typename INDEX_TYPE,
          int RESIDUAL_FORM >
class EquilibriumReactions
{
public:
  using RealType = REAL_TYPE;
  using IntType = INT_TYPE;
  using IndexType = INDEX_TYPE;

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
  enforceEquilibrium( RealType const & temperature,
                      PARAMS_DATA const & params,
                      ARRAY_1D_TO_CONST const & speciesConcentration0,
                      ARRAY_1D & speciesConcentration );

};


} // namespace bulkGeneric
} // namespace hpcReact

#include "EquilibriumReactions_impl.hpp"