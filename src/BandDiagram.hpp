#pragma once

// =============================================================================
// BandDiagram.hpp
// -----------------------------------------------------------------------------
// Renders a dynamic energy-band diagram:
//   * Valence band edge  E_v
//   * Conduction band edge  E_c
//   * Fermi level  E_f   (dashed line)
//   * Donor level (below E_c)  for n-type
//   * Acceptor level (above E_v)  for p-type
//   * Fermi-Dirac occupation f(E) as a curve overlaid on the right side
// =============================================================================

#include <SFML/Graphics.hpp>

#include "PhysicsEngine.hpp"


class BandDiagram {
public:
    BandDiagram(const sf::Vector2f& topLeft, const sf::Vector2f& size);

    // Renders the diagram.  `font` is required for axis + legend labels.
    void draw(sf::RenderTarget&   target,
              const PhysicsEngine& physics,
              const sf::Font&      font) const;

private:
    // Convert an energy value [eV] to a screen Y coordinate.
    float energyToY(double E) const;

    sf::Vector2f m_topLeft;
    sf::Vector2f m_size;

    // Energy axis range, in eV relative to E_v.
    float m_Emin = -0.30f;
    float m_Emax =  1.60f;
};
