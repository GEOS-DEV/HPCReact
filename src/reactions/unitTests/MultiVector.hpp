#pragma once

#include <vector>
#include <stdexcept>

template <typename T, size_t N>
class MultiVector {
public:

  // Default constructor: initializes dimensions to zero
  MultiVector() : m_size(0) {
    for (size_t i = 0; i < N; ++i) {
      m_dimensions[i] = 0;
    }
  }

  // Constructor: takes dimensions as variadic arguments
  template <typename... Dims>
  MultiVector(Dims... dims)
    : m_size(dims...) {
    static_assert(sizeof...(Dims) == N, "Incorrect number of dimensions.");
    size_t dimensions[N] = {static_cast<size_t>(dims)...};
    for (size_t i = 0; i < N; ++i) {
      m_dimensions[i] = dimensions[i];
    }
    m_data.resize(m_size);
  }

  // Constructor: takes dimensions as a C-style array
  MultiVector(const size_t (&dimensions)[N])
    : m_size(computeSize(dimensions)) {
    for (size_t i = 0; i < N; ++i) {
      m_dimensions[i] = dimensions[i];
    }
    m_data.resize(m_size);
  }

MultiVector(std::initializer_list<T> list)
    : m_size(list.size()), m_data(list) {
    static_assert(N == 1, "Initializer list constructor is only valid for 1D arrays.");
    m_dimensions[0] = list.size();
}

  // Resize function: allows changing dimensions and resizes underlying storage
  void resize(const size_t (&dimensions)[N]) {
    for (size_t i = 0; i < N; ++i) {
      m_dimensions[i] = dimensions[i];
    }
    m_size = computeSize(dimensions);
    m_data.resize(m_size);
  }

  // Operator() for element access
  template <typename... Indices>
  T& operator()(Indices... indices) {
    static_assert(sizeof...(Indices) == N, "Incorrect number of indices.");
    size_t index = computeIndex({static_cast<size_t>(indices)...});
    return m_data[index];
  }

  // Const version of operator()
  template <typename... Indices>
  const T& operator()(Indices... indices) const {
    static_assert(sizeof...(Indices) == N, "Incorrect number of indices.");
    size_t index = computeIndex({static_cast<size_t>(indices)...});
    return m_data[index];
  }

  // Operator() for element access
  template <typename... Indices>
  T& operator[](Indices... indices) {
    static_assert(sizeof...(Indices) == N, "Incorrect number of indices.");
    size_t index = computeIndex({static_cast<size_t>(indices)...});
    return m_data[index];
  }

  // Const version of operator()
  template <typename... Indices>
  const T& operator[](Indices... indices) const {
    static_assert(sizeof...(Indices) == N, "Incorrect number of indices.");
    size_t index = computeIndex({static_cast<size_t>(indices)...});
    return m_data[index];
  }

private:
  // Helper function to compute the total size from dimensions
  size_t computeSize(const size_t (&dimensions)[N]) const {
    size_t size = 1;
    for (size_t dim : dimensions) {
      size *= dim;
    }
    return size;
  }

  // Helper function to compute 1D index from multidimensional indices
  size_t computeIndex(const std::array<size_t, N>& indices) const {
    size_t index = 0;
    size_t multiplier = 1;
    for (size_t i = N; i-- > 0;) {
      if (indices[i] >= m_dimensions[i]) {
        throw std::out_of_range("Index out of range.");
      }
      index += indices[i] * multiplier;
      multiplier *= m_dimensions[i];
    }
    return index;
  }

  size_t m_dimensions[N];
  size_t m_size;
  std::vector<T> m_data;
};

template < int NDIM >
using RealMultiVector = MultiVector< double, NDIM >;

template < int NDIM >
using RealConstMultiVector = MultiVector< double const, NDIM >;
