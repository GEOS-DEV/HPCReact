#pragma once

#include "macros.hpp"
#include <cmath>

/**
 * @file logMath.hpp
 * @brief Base-selectable logarithm and exponential helpers with optional fast CUDA intrinsics.
 *
 * @details
 * This header provides a small, header-only utility for computing `log` and `exp` in either:
 * - natural base (e), or
 * - base-10 (ten),
 * and optionally using CUDA "fast math" intrinsics on device code.
 *
 * Key points:
 * - For `float` on CUDA device (`__CUDA_ARCH__`), `MathMode::fast` selects intrinsics such as
 *   `__logf`, `__log10f`, `__expf`, `__exp10f` when available.
 * - For host code (and for non-fast mode), it falls back to standard `<cmath>` functions
 *   like `::log`, `::log10`, `::exp`, `::expf`, `::logf`, `::log10f`.
 * - Some libcs (glibc) provide `exp10/exp10f` as extensions. When unavailable, base-10
 *   exponentiation is implemented via `exp(x * ln(10))`.
 *
 * @note
 * The functions here assume their usual mathematical domains:
 * - `log(x)` requires `x > 0`.
 * - `exp(x)` is defined for all real `x` but may overflow for large `x`.
 *
 * @warning
 * When `MathMode::fast` is used on CUDA device, results may differ slightly from the
 * accurate host/standard-library implementations due to reduced precision/approximations.
 */

namespace hpcReact
{

/**
 * @def HPCREACT_HAS_EXP10
 * @brief Indicates whether `exp10/exp10f` are available as libc extensions.
 *
 * @details
 * glibc commonly provides `exp10` and `exp10f` as non-standard extensions. This macro is used to
 * select a direct call to those functions on host when available; otherwise a mathematically
 * equivalent fallback is used:
 * \f[
 *   10^x = e^{x \ln(10)}
 * \f]
 */
#if defined(__GLIBC__)  // glibc provides exp10/exp10f as an extension
  #define HPCREACT_HAS_EXP10 1
#else
  #define HPCREACT_HAS_EXP10 0
#endif

/**
 * @brief Enumerates the logarithm base used by the helper.
 */
enum class LogBase
{
  e,   /**< Natural logarithm base (ln). */
  ten  /**< Base-10 logarithm. */
};

/**
 * @brief Enumerates the math evaluation mode.
 *
 * @details
 * `fast` is only used for `float` on CUDA device code. On host (or for `double`), the
 * standard library functions are used.
 */
enum class MathMode
{
  accurate, /**< Prefer standard-library implementations. */
  fast      /**< Prefer CUDA fast intrinsics for `float` on device when available. */
};

/**
 * @brief Base-selectable logarithm and exponential functions.
 *
 * @tparam Base Selects the logarithm/exp base:
 *   - `LogBase::e` for natural base
 *   - `LogBase::ten` for base-10
 * @tparam Mode Selects evaluation mode. See @ref MathMode.
 *
 * @details
 * Provides:
 * - `log(x)` returning \f$\ln(x)\f$ or \f$\log_{10}(x)\f$ depending on `Base`
 * - `exp(x)` returning \f$e^x\f$ or \f$10^x\f$ depending on `Base`
 *
 * The overload set covers `float` and `double`.
 */
template< LogBase Base, MathMode Mode = MathMode::accurate >
struct LogExp
{
  /**
   * @brief Compile-time constant for \f$\ln(10)\f$.
   *
   * @tparam T Floating-point type (typically `float` or `double`).
   * @return \f$\ln(10)\f$ as type `T`.
   *
   * @details
   * The literal is provided with long-double precision and cast to `T`.
   */
  template< class T >
  HPCREACT_HOST_DEVICE
  static constexpr T ln10()
  {
    return T(2.3025850929940456840179914546843642L);
  }

  /**
   * @brief Conversion factor for derivatives with respect to `log(x)` vs `ln(x)`.
   *
   * @tparam T Floating-point type (typically `float` or `double`).
   * @return A constant factor depending on @p Base.
   *
   * @details
   * This returns the scalar \f$c\f$ such that:
   * \f[
   *   \frac{d}{d(\log_b x)} = c \, \frac{d}{d(\ln x)}
   * \f]
   * where:
   * - if `Base == e`, then \f$c = 1\f$
   * - if `Base == ten`, then \f$c = \ln(10)\f$
   *
   * Equivalently, since \f$\log_{10}(x) = \ln(x) / \ln(10)\f$:
   * \f[
   *   \frac{d}{d(\log_{10} x)} = \ln(10) \, \frac{d}{d(\ln x)}
   * \f]
   */
  template< typename T >
  HPCREACT_HOST_DEVICE
  static constexpr inline T dWrtLogConst()
  {
    if constexpr ( Base == LogBase::e ) { return T(1.0); }
    else                                { return T(ln10<T>()); }
  }


  /**
   * @brief Adjusts a derivative value to account for logarithm base.
   *
   * @tparam T Floating-point type (typically `float` or `double`).
   * @param[in,out] x The derivative value to adjust.
   *
   * @details
   * This function multiplies the input `x` by the appropriate factor so that
   * it represents a derivative with respect to `log_{10}(x)` instead of `ln(x)`.
   * - If `Base == e`, `x` is unchanged.
   * - If `Base == ten`, `x` is multiplied by \f$\ln(10)\f$.
   */
  template< typename T >
  HPCREACT_HOST_DEVICE
  static constexpr inline void dWrtLogScale( T & x )
  {
    if constexpr ( Base == LogBase::ten ) { x = x * ln10<T>(); }
  }

  // ***** log function *********************************************************************************************

  /**
   * @brief Logarithm in the selected base for `double`.
   *
   * @param x Input value (must be > 0).
   * @return `ln(x)` if `Base == LogBase::e`, otherwise `log10(x)`.
   */
  HPCREACT_HOST_DEVICE
  static constexpr inline double log( double const x )
  {
    if constexpr ( Base == LogBase::e ) { return ::log( x ); }
    else                                { return ::log10( x ); }
  }

  /**
   * @brief Logarithm in the selected base for `float`.
   *
   * @param x Input value (must be > 0).
   * @return `ln(x)` if `Base == LogBase::e`, otherwise `log10(x)`.
   *
   * @details
   * On CUDA device code (`__CUDA_ARCH__`) and when `Mode == MathMode::fast`,
   * uses CUDA intrinsics:
   * - `__logf` / `__log10f`
   * Otherwise falls back to `<cmath>`:
   * - `::logf` / `::log10f`
   */
  HPCREACT_HOST_DEVICE
  static inline float log( float const x )
  {
#if defined(__CUDA_ARCH__)
    if constexpr ( Mode == MathMode::fast )
    {
      if constexpr ( Base == LogBase::e ) { return __logf( x ); }
      else                                { return __log10f( x ); }
    }
#endif
    if constexpr ( Base == LogBase::e ) { return ::logf( x ); }
    else                                { return ::log10f( x ); }
  }

  // ***** exp function *********************************************************************************************

  /**
   * @brief Exponential in the selected base for `double`.
   *
   * @param x Exponent.
   * @return `exp(x)` if `Base == LogBase::e`, otherwise `10^x`.
   *
   * @details
   * For base-10:
   * - If `HPCREACT_HAS_EXP10` is true, uses `::exp10(x)` (glibc extension).
   * - Otherwise uses `::exp(x * ln(10))`.
   */
  HPCREACT_HOST_DEVICE
  static inline double exp( double const x )
  {
    if constexpr ( Base == LogBase::e )
    {
      return ::exp( x );
    }
    else
    {
#if HPCREACT_HAS_EXP10
      return ::exp10( x );
#else
      return ::exp( x * ln10<double>() );
#endif
    }
  }

  /**
   * @brief Exponential in the selected base for `float`.
   *
   * @param x Exponent.
   * @return `exp(x)` if `Base == LogBase::e`, otherwise `10^x`.
   *
   * @details
   * On CUDA device code (`__CUDA_ARCH__`) and when `Mode == MathMode::fast`,
   * uses CUDA intrinsics:
   * - `__expf` / `__exp10f`
   *
   * Otherwise:
   * - For base-e uses `::expf(x)`
   * - For base-10 uses `::exp10f(x)` if available (glibc extension), else `::expf(x * ln(10))`.
   */
  HPCREACT_HOST_DEVICE
  static inline float exp( float const x )
  {
#if defined(__CUDA_ARCH__)
    if constexpr ( Mode == MathMode::fast )
    {
      if constexpr ( Base == LogBase::e ) { return __expf( x ); }
      else                                { return __exp10f( x ); }
    }
#endif
    if constexpr ( Base == LogBase::e )
    {
      return ::expf( x );
    }
    else
    {
#if HPCREACT_HAS_EXP10
      return ::exp10f( x );
#else
      return ::expf( x * ln10<float>() );
#endif
    }
  }
}; // struct LogExp

/**
 * @def HPC_REACT_LOG_TYPE
 * @brief Compile-time selection of the default `logmath` alias.
 *
 * @details
 * Values:
 * - `0`: natural log/exp (base-e), fast mode
 * - `1`: base-10 log/exp, fast mode (default)
 *
 * This macro is intended to allow build-time selection of the logarithm base used in code
 * that relies on the `hpcReact::logmath` alias.
 */
#if !defined(HPC_REACT_LOG_TYPE)
#define HPC_REACT_LOG_TYPE 1
#endif

/**
 * @brief Default log/exp helper type selected by `HPC_REACT_LOG_TYPE`.
 *
 * @details
 * - `HPC_REACT_LOG_TYPE == 0` selects `LogBase::e`
 * - `HPC_REACT_LOG_TYPE == 1` selects `LogBase::ten`
 *
 * Both selections currently use `MathMode::fast`. Note that `MathMode::fast` only changes
 * behavior for `float` on CUDA device code.
 */
#if HPC_REACT_LOG_TYPE == 0
  using logmath = LogExp< LogBase::e, MathMode::fast >;
#elif HPC_REACT_LOG_TYPE == 1
  using logmath = LogExp< LogBase::ten, MathMode::fast >;
#endif

}  // namespace hpcReact
