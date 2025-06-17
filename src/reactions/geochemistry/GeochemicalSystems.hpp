#pragma once

#include "Carbonate.hpp"
#include "Ultramafics.hpp"

namespace hpcReact
{

namespace geochemistry
{
  using systemTypes = std::variant< ultramaficSystemType,
                                    carbonateSystemType,
                                    carbonateSystemAllEquilibriumType >;

} // namespace geochemistry
} // namespace hpcReact
