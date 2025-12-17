
#pragma once

#include "macros.hpp"
#include <cmath>



namespace hpcReact
{

#if defined(__GLIBC__)  // glibc provides exp10/exp10f as an extension
  #define HPCREACT_HAS_EXP10 1
#else
  #define HPCREACT_HAS_EXP10 0
#endif

enum class LogBase { e, ten };
enum class MathMode { accurate, fast };

template< LogBase Base, MathMode Mode = MathMode::accurate >
struct LogExp
{

  template<class T>
  HPCREACT_HOST_DEVICE static constexpr 
  T ln10()
  {
    return T(2.3025850929940456840179914546843642L);
  }

  template< typename T >
  HPCREACT_HOST_DEVICE
  static constexpr inline 
  T dLogConst()
  {
    if constexpr ( Base == LogBase::e ) { return T(1.0); }
    else                                { return T(1.0/ln10()); }
  }



// ***** log function *********************************************************************************************
  HPCREACT_HOST_DEVICE 
  static inline 
  double log( double const x )
  {
    if constexpr ( Base == LogBase::e ) { return ::log(x); }
    else                                { return ::log10(x); }
  }
  
  HPCREACT_HOST_DEVICE 
  static inline 
  float log( float const x )
  {
#if defined(__CUDA_ARCH__)
    if constexpr ( Mode == MathMode::fast )
    {
      if constexpr ( Base == LogBase::e) { return __logf(x); }
      else                               { return __log10f(x); }
    }
#endif
    if constexpr ( Base == LogBase::e ) { return ::logf(x); }
    else                                { return ::log10f(x); }
  }

// ***** exp function *********************************************************************************************
  HPCREACT_HOST_DEVICE 
  static inline 
  double exp( double const x )
  {
    if constexpr ( Base == LogBase::e) 
    { 
      return ::exp(x); 
    }
    else
    { 
      #if HPCREACT_HAS_EXP10
        return ::exp10(x);
      #else
        return ::exp( x * ln10<double>() ); ; 
      #endif
    }
  }


  HPCREACT_HOST_DEVICE 
  static inline 
  float exp( float const x )
  {
#if defined(__CUDA_ARCH__)
    if constexpr ( Mode == MathMode::fast )
    {
      if constexpr ( Base == LogBase::e ) { return __expf(x); }
      else                                { return __exp10f(x); }
    }
#endif
    if constexpr ( Base == LogBase::e) 
    { 
      return ::expf(x); 
    }
    else
    { 
#if HPCREACT_HAS_EXP10
      return ::exp10f(x);
#else
      return ::expf( x * ln10<float>() ); 
#endif
    }
  }
}; // struct LogExp


#if !defined(HPC_REACT_LOG_TYPE)
#define HPC_REACT_LOG_TYPE 1
#endif

#if HPC_REACT_LOG_TYPE == 0
  using logmath = LogExp< LogBase::e, MathMode::fast >;
#elif HPC_REACT_LOG_TYPE == 1
  using logmath = LogExp< LogBase::ten, MathMode::fast >;
#endif

}  // namespace hpcReact