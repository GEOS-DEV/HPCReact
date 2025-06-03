#pragma once

#include "macros.hpp"
#include "DirectSystemSolve.hpp"
#include <cmath>

namespace hpcReact 
{

namespace nonlinearSolvers
{

namespace internal
{

template< int N >
HPCREACT_HOST_DEVICE 
double norm( double const (&r)[N]) 
{
    double sum = 0.0;
    for ( int i = 0; i < N; ++i ) sum += r[i] * r[i];
    return ::sqrt( sum );
}


template< int N >
HPCREACT_HOST_DEVICE 
void add(double (&x)[N], double const (&dx)[N]) 
{
    for (int i = 0; i < N; ++i) x[i] += dx[i];
}

template< int N >
HPCREACT_HOST_DEVICE 
void scale( double (&x)[N], double const value ) 
{
    for (int i = 0; i < N; ++i) x[i] *= value;
}

}

namespace utils
{   
template< int N >
HPCREACT_HOST_DEVICE
void print(double const (&J)[N][N], double const (&r)[N], double const (&dx)[N])
{
    printf("=======================\nJacobian matrix:\n=======================\n");
    printf("     RowID           ColID                       Value\n");    for ( int i = 0; i < N; ++i )
    {
        for ( int j = 0; j < N; ++j )
        {
            printf("%10d%16d%27.16e\n", i, j, J[i][j]);
        }
    }

    printf("\n=======================\nSystem right-hand side:\n=======================\n");
    printf("     RowID                       Value\n");

    for ( int i = 0; i < N; ++i )
    {
        printf("%10d%27.16e\n", i, r[i]);
    }
    
    printf("\n=======================\nDelta update vector:\n=======================\n");
    printf("     RowID                       Value\n");
    
    for ( int i = 0; i < N; ++i )
    {
        printf("%10d%27.16e\n", i, dx[i]);
    }
}
}

template< int N,
          typename REAL_TYPE, 
          typename ResidualFunc, 
          typename JacobianFunc >
HPCREACT_HOST_DEVICE 
void newtonRaphson( REAL_TYPE (&x)[N], 
                    ResidualFunc computeResidual,
                    JacobianFunc computeJacobian,
                    int maxIters = 25,
                    double tol = 1e-10 )
{
    REAL_TYPE residual[N];
    REAL_TYPE dx[N];
    REAL_TYPE jacobian[N][N];

    for ( int iter = 0; iter < maxIters; ++iter ) 
    {
        computeResidual( x, residual );

        if ( internal::norm< N >( residual ) < tol ) return;

        computeJacobian( x, jacobian );

        solveNxN_pivoted< REAL_TYPE, N >( jacobian, residual, dx );

        internal::add< N >( x, dx );
    }
}

template< int N,
          typename REAL_TYPE, 
          typename FUNCTION_TYPE > 
HPCREACT_HOST_DEVICE 
bool newtonRaphson( REAL_TYPE (&x)[N], 
                    FUNCTION_TYPE computeResidualAndJacobian,
                    int maxIters = 12,
                    double tol = 1e-10,
                    bool const do_print = false )
{
    REAL_TYPE residual[N]{};
    REAL_TYPE dx[N]{};
    REAL_TYPE jacobian[N][N]{};

     for ( int iter = 0; iter < maxIters; ++iter ) 
    {
        computeResidualAndJacobian( x, residual, jacobian );

        double const norm = internal::norm< N >( residual );

        printf( "--Iter %d: Residual norm = %.12e\n", iter, norm );

        if ( norm < tol ) {
            printf( "--Converged.\n" );
            return true;
        }
        internal::scale<N>( residual, -1.0);

        if( do_print )
        {
            utils::print( jacobian, residual, dx );
        }
       
        solveNxN_pivoted< REAL_TYPE, N >( jacobian, residual, dx );
        internal::add< N >( x, dx );
        
    }
    
    printf( "--Newton solver error: Max iterations reached without convergence.\n" );

    return false;
}

}
}