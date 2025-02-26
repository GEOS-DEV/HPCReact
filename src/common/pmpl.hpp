#pragma once

#include "common/macros.hpp"

#include <utility>


namespace hpcReact
{
#if defined(HPCREACT_USE_DEVICE)
  #if defined(HPCREACT_USE_CUDA)
#define deviceMalloc( PTR, BYTES ) cudaMalloc( PTR, BYTES );
#define deviceMallocManaged( PTR, BYTES ) cudaMallocManaged( PTR, BYTES );
#define deviceDeviceSynchronize() cudaDeviceSynchronize();
#define deviceMemCpy( DST, SRC, BYTES, KIND ) cudaMemcpy( DST, SRC, BYTES, KIND );
#define deviceFree( PTR ) cudaFree( PTR );
  #elif defined(HPCREACT_USE_HIP)
#define deviceMalloc( PTR, BYTES ) hipMalloc( PTR, BYTES );
#define deviceMallocManaged( PTR, BYTES ) hipMallocManaged( PTR, BYTES );
#define deviceDeviceSynchronize() hipDeviceSynchronize();
#define deviceMemCpy( DST, SRC, BYTES, KIND ) hipMemcpy( DST, SRC, BYTES, KIND );
#define deviceFree( PTR ) hipFree( PTR );
  #endif
#endif

/**
 * @namespace shiva::pmpl
 * @brief The pmpl namespace contains all of the pmpl classes and functions
 * used to provide a portablity layer in unit testing.
 */
namespace pmpl
{

/**
 * @brief This function checks if two floating point numbers are equal within a
 * tolerance.
 * @tparam REAL_TYPE This is the type of the floating point numbers to compare.
 * @param a This is the first floating point number to compare.
 * @param b This is the second floating point number to compare.
 * @param tolerance This is the tolerance to use when comparing the two numbers.
 * @return This returns true if the two numbers are equal within the tolerance.
 */
template< typename REAL_TYPE >
static constexpr bool check( REAL_TYPE const a, REAL_TYPE const b, REAL_TYPE const tolerance )
{
  return ( a - b ) * ( a - b ) < tolerance * tolerance;
}


/**
 * @brief This function provides a generic kernel execution mechanism that can
 * be called on either host or device.
 * @tparam LAMBDA The type of the lambda function to execute.
 * @param func The lambda function to execute.
 */
template< typename LAMBDA >
HPCREACT_GLOBAL void genericKernel( LAMBDA func )
{
  func();
}

/**
 * @brief This function provides a wrapper to the genericKernel function.
 * @tparam LAMBDA The type of the lambda function to execute.
 * @param func The lambda function to execute.
 *
 * This function will execute the lambda through a kernel launch of
 * genericKernel.
 */
template< typename LAMBDA >
void genericKernelWrapper( LAMBDA && func )
{
#if defined(HPCREACT_USE_DEVICE)
  genericKernel << < 1, 1 >> > ( std::forward< LAMBDA >( func ) );
#else
  genericKernel( std::forward< LAMBDA >( func ) );
#endif
}



/**
 * @brief This function provides a generic kernel execution mechanism that can
 * be called on either host or device.
 * @tparam DATA_TYPE The type of the data pointer.
 * @tparam LAMBDA The type of the lambda function to execute.
 * @param func The lambda function to execute.
 * @param data A general data pointer to pass to the lambda function that should
 * hold all data required to execute the lambda function, aside from what is
 * captured.
 */
template< typename DATA_TYPE, typename LAMBDA >
HPCREACT_GLOBAL void genericKernel( LAMBDA func, DATA_TYPE * const data )
{
  func( data );
}

/**
 * @brief This function provides a wrapper to the genericKernel function.
 * @tparam DATA_TYPE The type of the data pointer.
 * @tparam LAMBDA The type of the lambda function to execute.
 * @param N The size of the data array.
 * @param hostData The data pointer to pass to the lambda function.
 * @param func The lambda function to execute.
 *
 * This function will allocate the data pointer on the device, execute the
 * lambda through a kernel launch of genericKernel, and then synchronize the
 * device.
 */
template< typename DATA_TYPE, typename LAMBDA >
void genericKernelWrapper( int const N, DATA_TYPE * const hostData, LAMBDA && func )
{

#if defined(HPCREACT_USE_DEVICE)
  DATA_TYPE * deviceData;
  deviceMalloc( &deviceData, N * sizeof(DATA_TYPE) );
  deviceMemCpy( deviceData, hostData, N * sizeof(DATA_TYPE), cudaMemcpyHostToDevice );
  genericKernel <<< 1, 1 >>> ( std::forward< LAMBDA >( func ), deviceData );
  deviceDeviceSynchronize();
  deviceMemCpy( hostData, deviceData, N * sizeof(DATA_TYPE), cudaMemcpyDeviceToHost );
  deviceFree( deviceData );
#else
  HPCREACT_UNUSED_VAR( N );
  genericKernel( std::forward< LAMBDA >( func ), hostData );
#endif
}

/**
 * @brief convenience function for allocating data allocated on a pointer
 * @tparam DATA_TYPE The type of the data pointer.
 * @param data The data pointer to deallocate.
 */
template< typename DATA_TYPE >
HPCREACT_HOST_DEVICE void deallocateData( DATA_TYPE * data )
{
#if defined(HPCREACT_USE_DEVICE)
  deviceFree( data );
#else
  delete[] data;
#endif
}

} // namespace pmpl
} // namespace hpcReact
