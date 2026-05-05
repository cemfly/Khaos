#pragma once

// =============================================================================
// CrystalView.hpp -- Phase 5 + 5b + 6: textured 2D crystal panel
// -----------------------------------------------------------------------------
// Renders the crystal lattice + free carriers + (carrier OR thermal) heatmap
// + Lorentz vector field + (BJT) region overlay into an offscreen
// sf::RenderTexture so the result can be embedded in an ImGui dockable
// window.
//
// Layered composition (back -> front):
//
//   1. Background
//   2. Heatmap (carrier OR thermal)         -- toggleable
//   3. BJT region tint (E / B / C)          -- only when DeviceMode=NpnBjt
//   4. Lattice atoms + bonds                 -- always
//   5. Free carriers                         -- always
//   6. Vector field arrows (Lorentz)         -- toggleable
// =============================================================================

#include <SFML/Graphics.hpp>

#include <random>
#include <vector>

#include "DriftDiffusion.hpp"
#include "PhysicsEngine.hpp"


// -----------------------------------------------------------------------------
// Heatmap mode (None / Carriers / Thermal). The user picks one in the UI.
// -----------------------------------------------------------------------------
enum class HeatmapMode : int {
    None     = 0,
    Carriers = 1,
    Thermal  = 2,
};


class CrystalView {
public:
    explicit CrystalView(unsigned int textureSize = 720);

    // Sync atoms / carriers from the current physics state.
    void rebuild(const PhysicsEngine& physics);

    // Per-frame carrier random walk + Lorentz rotation.
    void update(float dt, const PhysicsEngine& physics);

    // Composite the layered scene to the internal RenderTexture.
    void render(const PhysicsEngine&  physics,
                const DriftDiffusion& dd,
                HeatmapMode           heatmapMode,
                bool                  showVectorField);

    [[nodiscard]] const sf::RenderTexture& renderTexture() const noexcept {
        return m_rt;
    }
    [[nodiscard]] unsigned int size() const noexcept {
        return m_rt.getSize().x;
    }

private:
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

    void resampleCarriers(const PhysicsEngine& physics);
    [[nodiscard]] int targetCarrierCount(double concentration) const noexcept;

    void drawCarrierHeatmap(const DriftDiffusion& dd);
    void drawThermalHeatmap(const DriftDiffusion& dd);
    void drawBjtRegions    (const DriftDiffusion& dd);
    void drawLattice       (const PhysicsEngine& physics);
    void drawCarriers      ();
    void drawVectorField   (const PhysicsEngine& physics);

    // ---- State -------------------------------------------------------------
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
