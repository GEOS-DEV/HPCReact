
#ifndef HPCREACT_MACROS_HPP
#define HPCREACT_MACROS_HPP

#if defined( __CUDACC__ ) || defined( __HIPCC__ )
#define HPCREACT_USE_DEVICE
#define HPCREACT_HOST_DEVICE __host__ __device__
#else
#define HPCREACT_HOST_DEVICE
#endif


#if defined( HPCREACT_USE_DEVICE )
#define HPCREACT_GLOBAL __global__
#else
#define HPCREACT_GLOBAL
#endif

/// This macro is used to ignore warnings that that a variable is
/// unused.
#define HPCREACT_UNUSED_VAR( ... ) (void)( __VA_ARGS__ )

#endif // HPCREACT_MACROS_HPP