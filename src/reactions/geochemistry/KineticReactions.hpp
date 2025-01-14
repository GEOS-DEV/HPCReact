#pragma once

#include "ReactionsBase.hpp"
#include "common/macros.hpp"

namespace hpcReact
{
namespace geochemistry
{

template< typename REAL_TYPE,
          typename REAL_DATA_ARRAY_1D_VIEW_TYPE,
          typename REAL_CONST_DATA_ARRAY_1D_VIEW_TYPE,
          typename INT_TYPE,
          typename INDEX_TYPE >
class KineticReactions : public ReactionsBase< REAL_TYPE,
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


  static constexpr RealType RConst = 0; //constants::gasConstant;

  template< typename PARAMS_DATA >
  static HPCREACT_HOST_DEVICE inline void
  computeReactionRates( RealType const & temperature,
                        PARAMS_DATA const & params,
                        RealConstDataArrayView1d & primarySpeciesConcentration,
                        RealConstDataArrayView1d & secondarySpeciesConcentration,
                        RealDataArrayView1d & reactionRates );


};

} // namespace geochemistry
} // namespace hpcReact

#include "common/macrosCleanup.hpp"
