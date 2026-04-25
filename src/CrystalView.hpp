#pragma once

// =============================================================================
// CrystalView.hpp -- Phase 5 + 5b: textured 2D crystal panel
// -----------------------------------------------------------------------------
// Renders the 2D crystal lattice + free carriers + heatmap + Lorentz vector
// field into an offscreen sf::RenderTexture so the result can be embedded
// in an ImGui dockable window via ImGui::Image(const sf::RenderTexture&, ...).
//
// Layered composition (back -> front):
//
//   1. Background
//   2. Heatmap   (n(x,y) from DriftDiffusion)        -- toggleable
//   3. Lattice atoms + bonds                         -- always
//   4. Free carriers (random walk + Lorentz curl)    -- always
//   5. Vector field arrows (Lorentz deflection)      -- toggleable
//
// =============================================================================

#include <SFML/Graphics.hpp>

#include <random>
#include <vector>

#include "DriftDiffusion.hpp"
#include "PhysicsEngine.hpp"


class CrystalView {
public:
    explicit CrystalView(unsigned int textureSize = 720);

    // ---- State management --------------------------------------------------
    // Sync atom and carrier arrays with the current physics state. Should
    // be called whenever doping / material / temperature / optical state
    // changes meaningfully.
    void rebuild(const PhysicsEngine& physics);

    // Per-frame update of the carrier random walk + Lorentz rotation.
    void update(float dt, const PhysicsEngine& physics);

    // ---- Rendering ---------------------------------------------------------
    // Composites the layered scene into the internal RenderTexture.
    void render(const PhysicsEngine&  physics,
                const DriftDiffusion& dd,
                bool                  showHeatmap,
                bool                  showVectorField);

    // RenderTexture for ImGui::Image embedding. Stable reference -- the
    // texture is owned by this object.
    [[nodiscard]] const sf::RenderTexture& renderTexture() const noexcept {
        return m_rt;
    }

    // Texture size in pixels (square).
    [[nodiscard]] unsigned int size() const noexcept {
        return m_rt.getSize().x;
    }

private:
    // ---- Internal types ----------------------------------------------------
    enum class AtomType { Host, Donor, Acceptor };

    struct Atom {
        sf::Vector2f pos;
        AtomType     type;
    };

    struct Carrier {
        sf::Vector2f pos;
        sf::Vector2f vel;
        bool         electron;
        bool         optical;
    };

    // ---- Helpers -----------------------------------------------------------
    void resampleCarriers(const PhysicsEngine& physics);
    [[nodiscard]] int targetCarrierCount(double concentration) const noexcept;

    void drawHeatmap   (const DriftDiffusion& dd);
    void drawLattice   (const PhysicsEngine& physics);
    void drawCarriers  ();
    void drawVectorField(const PhysicsEngine& physics);

    // ---- Geometry / state --------------------------------------------------
    sf::RenderTexture m_rt;
    sf::Vector2f      m_topLeft;
    sf::Vector2f      m_size;
    int               m_rows;
    int               m_cols;
    float             m_cellW;
    float             m_cellH;

    std::vector<Atom>    m_atoms;
    std::vector<Carrier> m_carriers;

    mutable std::mt19937 m_rng;
};
