#pragma once

// =============================================================================
// Palette.hpp
// -----------------------------------------------------------------------------
// Shared visual utilities used by CrystalLattice, BandDiagram and UIPanel.
//
//   * Central colour palette -- a single source of truth for every hue shown
//     on screen, so the three views agree on what, say, a "phosphorus atom"
//     or a "thermal hole" looks like.
//   * Wavelength -> RGB mapping for the optical slider / arrow.
//   * makeText() factory that encapsulates the SFML 3 Text construction
//     (font is mandatory, no default constructor) and returns a configured
//     sf::Text in one line.
// =============================================================================

#include <SFML/Graphics.hpp>

#include <string>


namespace palette {

    // ---- Semantic colours --------------------------------------------------
    // Atom species
    inline const sf::Color Silicon   {120, 170, 220};
    inline const sf::Color Phosphorus{ 80, 220, 120};
    inline const sf::Color Boron     {230, 110, 110};

    // Thermal carriers
    inline const sf::Color Electron  {255, 230,  80};   // yellow
    inline const sf::Color Hole      {255,  80, 200};   // magenta

    // Photo-generated (optical) carriers
    inline const sf::Color ElectronOpt{140, 255, 255};  // cyan
    inline const sf::Color HoleOpt    {255, 180, 255};  // light pink

    // UI chrome
    inline const sf::Color PanelBg   { 28,  33,  45};
    inline const sf::Color PanelEdge { 80,  90, 110};
    inline const sf::Color ViewBg    { 18,  22,  30};
    inline const sf::Color WindowBg  { 10,  12,  18};
    inline const sf::Color TextLight {230, 235, 245};
    inline const sf::Color TextDim   {200, 205, 220};
    inline const sf::Color Accent    {255, 240, 150};
    inline const sf::Color FermiLine {255, 240,  90};

    // Band colours
    inline const sf::Color ConductionBand{230, 140, 100};
    inline const sf::Color ValenceBand   {120, 180, 230};
    inline const sf::Color DonorLevel    { 80, 230, 140};
    inline const sf::Color AcceptorLevel {240, 120, 120};

    // Transport readouts
    inline const sf::Color Conductivity {140, 255, 180};
    inline const sf::Color Hall         {255, 180, 120};
    inline const sf::Color Optical      {140, 230, 255};


    // -----------------------------------------------------------------------
    // Wavelength to RGB (pedagogical linear mapping across the visible
    // spectrum). UV (< 380 nm) and IR (> 780 nm) fall back to neutral grey.
    // -----------------------------------------------------------------------
    [[nodiscard]] inline sf::Color wavelengthColor(double lambda_nm) noexcept {
        if (lambda_nm < 380.0 || lambda_nm > 780.0) return {160, 160, 180};
        if (lambda_nm < 440.0)                      return {130, 100, 255};
        if (lambda_nm < 490.0)                      return { 90, 180, 255};
        if (lambda_nm < 510.0)                      return { 80, 255, 220};
        if (lambda_nm < 580.0)                      return {120, 255, 120};
        if (lambda_nm < 645.0)                      return {255, 220,  80};
        return                                             {255, 110,  90};
    }


    // -----------------------------------------------------------------------
    // sf::Text factory (SFML 3 requires a Font in the constructor).
    // -----------------------------------------------------------------------
    [[nodiscard]] inline sf::Text makeText(const sf::Font&     font,
                                           const std::string&  str,
                                           unsigned int        size,
                                           sf::Color           colour,
                                           const sf::Vector2f& pos)
    {
        sf::Text t(font, str, size);
        t.setFillColor(colour);
        t.setPosition(pos);
        return t;
    }

} // namespace palette
