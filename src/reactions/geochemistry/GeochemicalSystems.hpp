#pragma once

#include "Carbonate.hpp"
#include "Ultramafics.hpp"
#include "Forge.hpp"

#include <variant>

namespace hpcReact
{

namespace geochemistry
{
using systemTypes = std::variant< ultramaficSystemType,
                                  carbonateSystemType,
                                  carbonateSystemAllEquilibriumType,
                                  forgeSystemType >;

} // namespace geochemistry
} // namespace hpcReact
