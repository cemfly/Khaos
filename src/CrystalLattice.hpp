#pragma once

// =============================================================================
// CrystalLattice.hpp    (Phase 3)
// -----------------------------------------------------------------------------
// 2D schematic visualization of a Silicon crystal lattice.
//
// Phase 1 / 2 features
//   * Si atoms on a square lattice with bond lines
//   * Dopant atoms (P / B) replacing a fraction of Si sites
//   * Free electrons / holes performing a damped random walk
//
// Phase 3 additions
//   * Lorentz force  F = q (v x B)  applied to every carrier so that
//     electrons and holes curve in opposite directions when B != 0.
//   * Optical excess: when the physics engine reports Delta n > 0 (above-gap
//     illumination), additional electron-hole pairs are visualized.  The
//     dopant pattern is unchanged -- only the free-carrier cloud grows.
// =============================================================================

#include <SFML/Graphics.hpp>

#include <random>
#include <vector>

#include "PhysicsEngine.hpp"


class CrystalLattice {
public:
    CrystalLattice(const sf::Vector2f& topLeft,
                   const sf::Vector2f& size,
                   int rows,
                   int cols);

    // Regenerates dopant sites + carrier population for the current physics
    // state.  Call whenever doping type / concentration / optical state
    // changes.
    void rebuild(const PhysicsEngine& physics);

    // Advance the carrier random walk (and apply Lorentz rotation if B != 0).
    void update(float dt, const PhysicsEngine& physics);

    // Render the lattice + carriers.
    void draw(sf::RenderTarget& target) const;

private:
    enum class AtomType { Silicon, Phosphorus, Boron };

    struct Atom {
        sf::Vector2f pos;
        AtomType     type;
    };

    struct Carrier {
        sf::Vector2f pos;
        sf::Vector2f vel;
        bool         electron;   // true = electron, false = hole
        bool         optical;    // true => pair came from photon absorption
    };

    // Helpers
    void resampleCarriers(const PhysicsEngine& physics);
    int  targetCarrierCount(double concentration) const;

    // Geometry
    sf::Vector2f m_topLeft;
    sf::Vector2f m_size;
    int   m_rows;
    int   m_cols;
    float m_cellW;
    float m_cellH;

    // State
    std::vector<Atom>    m_atoms;
    std::vector<Carrier> m_carriers;

    mutable std::mt19937 m_rng;
};
