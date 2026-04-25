#pragma once

// =============================================================================
// BandView.hpp -- Phase 5: textured band diagram panel
// -----------------------------------------------------------------------------
// Renders the energy band diagram (E_v / E_c / E_f / dopant level) and the
// Fermi-Dirac distribution into an offscreen RenderTexture so it can be
// embedded inside an ImGui dockable window.
// =============================================================================

#include <SFML/Graphics.hpp>

#include "PhysicsEngine.hpp"


class BandView {
public:
    BandView(unsigned int width = 720, unsigned int height = 540);

    // Repaint the texture from the current physics state.
    void render(const PhysicsEngine& physics, const sf::Font& font);

    [[nodiscard]] const sf::RenderTexture& renderTexture() const noexcept {
        return m_rt;
    }

    [[nodiscard]] unsigned int width()  const noexcept { return m_rt.getSize().x; }
    [[nodiscard]] unsigned int height() const noexcept { return m_rt.getSize().y; }

private:
    [[nodiscard]] float energyToY(double E) const noexcept;

    sf::RenderTexture m_rt;
    float             m_Emin = -0.30f;
    float             m_Emax =  1.80f;
};
