#pragma once

#include "common/macros.hpp"
#include "common/CArrayWrapper.hpp"
#include "common/DirectSystemSolve.hpp"
#include <cmath>

#include<stdio.h>


namespace hpcReact
{
namespace bulkGeneric
{

template< typename REAL_TYPE,
          typename INT_TYPE,
          typename INDEX_TYPE,
          bool LOGE_CONCENTRATION,
          bool LOGE_RESIDUAL >
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
                                    ARRAY_2D & jacobian )
  {
    HPCREACT_UNUSED_VAR( temperature );
    constexpr int numSpecies = PARAMS_DATA::numSpecies;
    constexpr int numReactions = PARAMS_DATA::numReactions;
    
    // initialize the species concentration
    RealType speciesConcentration[numSpecies];
    for( IndexType i=0; i<numSpecies; ++i )
    {
      speciesConcentration[i] = speciesConcentration0[i];
      for( IndexType r=0; r<numReactions; ++r )
      {
        speciesConcentration[i] += params.stoichiometricMatrix( r, i ) * xi[r];
      }
    }
    
    // loop over reactions
    for( IndexType a=0; a<numReactions; ++a )
    {
      // get the equilibrium constant for this reaction
      RealType const Keq = params.equilibriumConstant( a );

      RealType forwardProduct = 1.0;
      RealType reverseProduct = 1.0;
      RealType dForwardProduct_dxi[numReactions] = {0.0};
      RealType dReverseProduct_dxi[numReactions] = {0.0};
      // loop over species
      for( IndexType i=0; i<numSpecies; ++i )
      {
        RealType const s_ai = params.stoichiometricMatrix( a, i );

        if( s_ai < 0.0 )
        {
          // forward reaction
          forwardProduct *= pow( speciesConcentration[i], -s_ai );
          // derivative of forward product with respect to xi
          for( IndexType b=0; b<numReactions; ++b )
          {
            dForwardProduct_dxi[b] += -s_ai / speciesConcentration[i] * params.stoichiometricMatrix( b, i );
          }
        }
        else if( s_ai > 0.0 )
        {
          // reverse reaction
          reverseProduct *= pow( speciesConcentration[i], s_ai );
          // derivative of reverse product with respect to xi
          for( IndexType b=0; b<numReactions; ++b )
          {
            dReverseProduct_dxi[b] += s_ai / speciesConcentration[i] * params.stoichiometricMatrix( b, i );
          }
        }
      }
      // compute the residual for this reaction
     residual[a] = forwardProduct - Keq * reverseProduct;
//     printf( "residual[%d] = %g - %g * %g = %g\n", a, forwardProduct, Keq, reverseProduct, residual[a] );

      // Finish the derivatives of the product terms with respect to xi
      for( IndexType b=0; b<numReactions; ++b )
      {
        dForwardProduct_dxi[b] *= forwardProduct;
        dReverseProduct_dxi[b] *= reverseProduct;
        // printf( "(%d) dForwardProduct_dxi[%d] = %f\n", a, b, dForwardProduct_dxi[b] );
        // printf( "(%d) dReverseProduct_dxi[%d] = %f\n", a, b, dReverseProduct_dxi[b] );
        // printf( "Keq = %f\n", Keq );
        jacobian( a, b ) = dForwardProduct_dxi[b] - Keq * dReverseProduct_dxi[b];
      }
    }
  }


  template< typename PARAMS_DATA,
  typename ARRAY_1D,
  typename ARRAY_1D_TO_CONST >
  static HPCREACT_HOST_DEVICE 
  void
  enforceEquilibrium( RealType const & temperature,
                                    PARAMS_DATA const & params,
                                    ARRAY_1D_TO_CONST const & speciesConcentration0,
                                    ARRAY_1D & speciesConcentration )
  {
    HPCREACT_UNUSED_VAR( temperature );
    constexpr int numSpecies = PARAMS_DATA::numSpecies;
    constexpr int numReactions = PARAMS_DATA::numReactions;
    double residual[numReactions] = { 0.0 };
    double xi[numReactions] = { 0.0 };
    double dxi[numReactions] = { 0.0 };
    CArrayWrapper< double, numReactions, numReactions > jacobian;
  
    REAL_TYPE residualNorm = 0.0;
    for( int k=0; k<10; ++k )
    {
      computeResidualAndJacobian( temperature,
                                  params,
                                  speciesConcentration0,
                                  xi,
                                  residual,
                                  jacobian );
      
      residualNorm = 0.0;
      for( int j = 0; j < numSpecies; ++j )
      {
        residualNorm += residual[j] * residual[j];
      }
      residualNorm = sqrt( residualNorm );
      if( residualNorm < 1.0e-14 )
      {
        break;
      }

      solveNxN_pivoted<double,numReactions>( jacobian.data, residual, dxi );

      // solve for the change in xi
      for( IndexType r=0; r<numReactions; ++r )
      {
        xi[r] -= dxi[r];
      }
    }

    for( IndexType i=0; i<numSpecies; ++i )
    {
      speciesConcentration[i] = speciesConcentration0[i];
      for( IndexType r=0; r<numReactions; ++r )
      {
        speciesConcentration[i] += params.stoichiometricMatrix( r, i ) * xi[r];
      }
    }
  }

};


} // namespace bulkGeneric
} // namespace hpcReact