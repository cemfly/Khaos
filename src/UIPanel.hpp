#pragma once

// =============================================================================
// UIPanel.hpp    (Phase 3)
// -----------------------------------------------------------------------------
// Interactive control panel for the simulator.
//
// Widgets
//   Row A (sliders)
//     * Temperature slider         (linear,        100 K .. 600 K)
//     * Doping concentration       (logarithmic,  1e14 .. 1e19 cm^-3)
//   Row B (sliders)
//     * Wavelength lambda          (linear,  300 .. 1500 nm)
//     * Magnetic field B           (linear,    0 ..   10 T)
//
//   Row C (buttons)
//     * Intrinsic / n-type / p-type
//     * Ionization: Full / Incomplete
//     * Light: On / Off
//     * Mobility: Matthiessen / Arora
//     * Export CSV
//
// Readouts (right-hand column)
//     T, E_g, E_f, E_photon, n_i, n, p, Delta n, mu_n, mu_p, sigma, rho,
//     R_H, ionization fraction.
// =============================================================================

#include <SFML/Graphics.hpp>

#include <memory>
#include <string>

#include "PhysicsEngine.hpp"


// -----------------------------------------------------------------------------
// Slider:  horizontal 1D handle with linear OR base-10-logarithmic mapping.
// -----------------------------------------------------------------------------
class Slider {
public:
    Slider(const sf::Vector2f& pos,
           float length,
           float minVal,
           float maxVal,
           float initVal,
           bool  logarithmic = false);

    void handleEvent(const sf::Event& ev, const sf::Vector2f& mousePos);

    void draw(sf::RenderTarget&  target,
              const std::string& label,
              const std::string& valueString,
              const sf::Font&    font) const;

    float getValue() const { return m_value; }
    bool  consumeChanged() { bool c = m_changed; m_changed = false; return c; }
    sf::FloatRect bounds() const;

private:
    float normalized() const;
    void  setFromNormalized(float t);

    sf::Vector2f m_pos;
    float m_length;
    float m_min;
    float m_max;
    float m_value;
    bool  m_log;
    bool  m_dragging = false;
    bool  m_changed  = false;
};


// -----------------------------------------------------------------------------
// Rectangular button with a latched "active" styling.
// -----------------------------------------------------------------------------
class Button {
public:
    Button(const sf::Vector2f& pos,
           const sf::Vector2f& size,
           std::string         label);

    void handleEvent(const sf::Event& ev, const sf::Vector2f& mousePos);
    void draw(sf::RenderTarget& target, const sf::Font& font) const;

    void setActive(bool a)          { m_active = a; }
    void setLabel (std::string s)   { m_label  = std::move(s); }
    bool isActive()    const        { return m_active;  }
    bool consumeClicked() { bool c = m_clicked; m_clicked = false; return c; }

private:
    sf::Vector2f m_pos;
    sf::Vector2f m_size;
    std::string  m_label;
    bool         m_active  = false;
    bool         m_clicked = false;
};


// -----------------------------------------------------------------------------
// UIPanel
// -----------------------------------------------------------------------------
class UIPanel {
public:
    UIPanel(const sf::Vector2f& pos, const sf::Vector2f& size);

    // Returns true when the visualization should rebuild.
    bool handleEvent(const sf::Event&     ev,
                     const sf::Vector2f&  mousePos,
                     PhysicsEngine&       physics);

    void draw(sf::RenderTarget&   target,
              const PhysicsEngine& physics,
              const sf::Font&      font) const;

private:
    void performCSVExport(const PhysicsEngine& physics);

    sf::Vector2f m_pos;
    sf::Vector2f m_size;

    // Sliders
    std::unique_ptr<Slider> m_tempSlider;
    std::unique_ptr<Slider> m_dopingSlider;
    std::unique_ptr<Slider> m_lambdaSlider;
    std::unique_ptr<Slider> m_bFieldSlider;

    // Doping type
    std::unique_ptr<Button> m_btnIntrinsic;
    std::unique_ptr<Button> m_btnNType;
    std::unique_ptr<Button> m_btnPType;

    // Mode toggles
    std::unique_ptr<Button> m_btnIonization;    // Full / Incomplete
    std::unique_ptr<Button> m_btnLight;         // On / Off
    std::unique_ptr<Button> m_btnMobility;      // Matthiessen / Arora

    // Action
    std::unique_ptr<Button> m_btnExportCSV;

    // Status message
    std::string       m_statusMessage;
    mutable sf::Clock m_statusClock;
    float             m_statusLifetime = 0.f;
};
