#pragma once

// =============================================================================
// BandView.hpp -- textured band diagram panel (flat + spatially bent modes)
//
//   Author : dex / cemfly-april2026
//   License: MIT
// -----------------------------------------------------------------------------
// Two render modes:
//
//   * Flat  : single-point band diagram derived from PhysicsEngine state
//             (E_v, E_c, E_f, donor / acceptor levels, Fermi-Dirac panel).
//
//   * Spatial : bends E_c(x), E_v(x), E_f along a horizontal cut taken
//               from the DriftDiffusion grid.  Drives the Phase 2 PN-
//               junction band-bending visualisation.
//
//                   E_c(x) = E_c0 - q psi(x)
//                   E_v(x) = E_v0 - q psi(x)
//
//               (Sze Eq. 2.10; Pierret Sec. 4.2.) The cut row is set by
//               the UI ("Cut y" slider). Below the bands we also overlay
//               the local |E|(x) curve, which is what feeds the
//               Chynoweth / Kane BTBT rate plots.
//
// The render texture is owned here so the view can be embedded inside an
// ImGui dockable window via ImGui::Image(view.renderTexture(), ...).
// =============================================================================

#include <SFML/Graphics.hpp>

#include <span>
#include <vector>

#include "DriftDiffusion.hpp"
#include "PhysicsEngine.hpp"


class BandView {
public:
    enum class Mode : std::uint8_t {
        Flat    = 0,    // legacy single-point diagram
        Spatial = 1,    // bent bands along a horizontal cut
    };

    BandView(unsigned int width = 720, unsigned int height = 540);

    // ---- Render: flat (legacy, single point) ----------------------------
    void render(const PhysicsEngine& physics, const sf::Font& font);

    // ---- Render: spatial cut (bent bands) -------------------------------
    //
    // Pulls Ec(x), Ev(x), |E|(x) along row j_cut from the drift-diffusion
    // grid and paints them onto the texture together with the constant
    // Fermi level. Uses *internal* pre-allocated buffers; this function
    // is on the hot path, so allocations happen only once (lazy resize
    // when the grid width changes).
    void renderSpatial(const PhysicsEngine& physics,
                       const DriftDiffusion& dd,
                       int                   j_cut,
                       const sf::Font&       font);

    void setMode(Mode m) noexcept { m_mode = m; }
    [[nodiscard]] Mode mode() const noexcept { return m_mode; }

    [[nodiscard]] const sf::RenderTexture& renderTexture() const noexcept {
        return m_rt;
    }
    [[nodiscard]] unsigned int width()  const noexcept { return m_rt.getSize().x; }
    [[nodiscard]] unsigned int height() const noexcept { return m_rt.getSize().y; }

private:
    [[nodiscard]] float energyToY(double E) const noexcept;
    [[nodiscard]] float xCellToPx(int    i, int W_cells,
                                  float xL, float xR) const noexcept;

    void ensureBuffersFor(int W_cells);

    sf::RenderTexture m_rt;
    Mode              m_mode = Mode::Spatial;

    float             m_Emin = -0.30f;
    float             m_Emax =  1.80f;

    // Lazy-allocated cut buffers (Phase 2 + 3). Grow only on first use
    // or when the dd grid width changes; never shrink.
    std::vector<float> m_psi_cut;
    std::vector<float> m_Ec_cut;
    std::vector<float> m_Ev_cut;
    std::vector<float> m_Efield_cut;
    std::vector<float> m_doping_cut;
    std::vector<float> m_phiN_cut;     // electron quasi-Fermi pot. (V)
    std::vector<float> m_phiP_cut;     // hole     quasi-Fermi pot. (V)
};
