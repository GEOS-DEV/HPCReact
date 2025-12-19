/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: (BSD-3-Clause)
 *
 * Copyright (c) 2025- Lawrence Livermore National Security LLC
 * All rights reserved
 *
 * See top level LICENSE files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

#pragma once



#if defined( __CUDACC__ )
#define HPCREACT_USE_DEVICE
#define HPCREACT_USE_CUDA
#elif defined( __HIPCC__ )
#define HPCREACT_USE_DEVICE
#define HPCREACT_USE_HIP
#endif


#if defined( HPCREACT_USE_DEVICE )
#define HPCREACT_GLOBAL __global__
#define HPCREACT_DEVICE __device__
#define HPCREACT_HOST_DEVICE __host__ __device__
#else
#define HPCREACT_GLOBAL
#define HPCREACT_DEVICE
#define HPCREACT_HOST_DEVICE
#endif

/// This macro is used to ignore warnings that that a variable is
/// unused.
#define HPCREACT_UNUSED_VAR( ... ) (void)( __VA_ARGS__ )


#if defined( __clang__ )
#define HPCREACT_NO_MISSING_BRACES( ... ) \
        _Pragma("clang diagnostic push") \
        _Pragma("clang diagnostic ignored \"-Wmissing-braces\"") \
        __VA_ARGS__ \
        _Pragma("clang diagnostic pop")

#define HPCREACT_NO_MISSING_BRACES_OPEN \
        _Pragma("clang diagnostic push") \
        _Pragma("clang diagnostic ignored \"-Wmissing-braces\"")
#define HPCREACT_NO_MISSING_BRACES_CLOSE \
        _Pragma("clang diagnostic pop")

#elif defined(__GNUC__)
#define HPCREACT_NO_MISSING_BRACES( ... ) \
        _Pragma("GCC diagnostic push") \
        _Pragma("GCC diagnostic ignored \"-Wmissing-braces\"") \
        __VA_ARGS__ \
        _Pragma("GCC diagnostic pop")

#define HPCREACT_NO_MISSING_BRACES_OPEN \
        _Pragma("GCC diagnostic push") \
        _Pragma("GCC diagnostic ignored \"-Wmissing-braces\"")
#define HPCREACT_NO_MISSING_BRACES_CLOSE \
        _Pragma("GCC diagnostic pop")

#elif defined(_MSC_VER)
#define HPCREACT_NO_MISSING_BRACES( ... ) \
        __pragma(warning(push)) \
        __pragma(warning(disable : 4351)) \
        __VA_ARGS__ \
        __pragma(warning(pop))

#define HPCREACT_NO_MISSING_BRACES_OPEN \
        __pragma(warning(push)) \
        __pragma(warning(disable : 4351))
#define HPCREACT_NO_MISSING_BRACES_CLOSE \
        __pragma(warning(pop))
#else
#define HPCREACT_NO_MISSING_BRACES( ... ) __VA_ARGS__ // No-op for unknown compilers
#define HPCREACT_NO_MISSING_BRACES_OPEN
#define HPCREACT_NO_MISSING_BRACES_CLOSE
#endif
