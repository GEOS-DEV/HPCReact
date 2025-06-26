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
