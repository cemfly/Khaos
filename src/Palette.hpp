#pragma once

// =============================================================================
// Palette.hpp -- shared visual constants
//
//   Author : dex / cemfly-april2026
//   License: MIT
// -----------------------------------------------------------------------------
// Single source of truth for every colour the simulator draws. Tuned for a
// "light scientific" aesthetic (COMSOL / MATLAB / Mathematica feel) -- white
// backgrounds, dark slate text, saturated but slightly desaturated accents
// so plots and overlays read clearly on light surfaces.
// =============================================================================

#include <SFML/Graphics.hpp>

#include <cstdint>
#include <string>


namespace palette {

    // ---- UI chrome (light scientific) -------------------------------------
    inline const sf::Color WindowBg   {232, 235, 240};   // outer SFML clear
    inline const sf::Color PanelBg    {248, 249, 252};   // ImGui panels
    inline const sf::Color ViewBg     {252, 252, 254};   // SFML render textures
    inline const sf::Color PanelEdge  {185, 195, 210};   // borders / frames
    inline const sf::Color TextLight  { 28,  34,  46};   // primary text
    inline const sf::Color TextDim    { 95, 105, 125};   // secondary text
    inline const sf::Color Accent     { 60, 110, 200};   // sliders / handles
    inline const sf::Color FermiLine  {184, 134,  11};   // dark goldenrod

    // ---- Atom species (saturated for visibility on white) -----------------
    inline const sf::Color Silicon    { 60, 110, 175};
    inline const sf::Color Phosphorus { 40, 130,  70};
    inline const sf::Color Boron      {180,  60,  60};

    // ---- Thermal carriers --------------------------------------------------
    inline const sf::Color Electron   {220, 130,  50};   // orange
    inline const sf::Color Hole       {150,  50, 180};   // deep purple

    // ---- Photo-generated carriers -----------------------------------------
    inline const sf::Color ElectronOpt{ 50, 130, 220};   // bright blue
    inline const sf::Color HoleOpt    {220, 100, 180};   // pink

    // ---- Band edges and dopant levels -------------------------------------
    inline const sf::Color ConductionBand{200,  85,  35};
    inline const sf::Color ValenceBand   { 30,  90, 170};
    inline const sf::Color DonorLevel    { 30, 140,  80};
    inline const sf::Color AcceptorLevel {200,  60,  60};

    // ---- Transport readouts ------------------------------------------------
    inline const sf::Color Conductivity  { 30, 140,  90};
    inline const sf::Color Hall          {200, 110,  40};
    inline const sf::Color Optical       { 50, 130, 220};


    // -----------------------------------------------------------------------
    // Wavelength -> RGB (visible spectrum); UV / IR fall back to mid grey.
    // -----------------------------------------------------------------------
    [[nodiscard]] inline sf::Color wavelengthColor(double lambda_nm) noexcept {
        if (lambda_nm < 380.0 || lambda_nm > 780.0) return {130, 140, 160};
        if (lambda_nm < 440.0)                      return {110,  60, 220};
        if (lambda_nm < 490.0)                      return { 50, 110, 220};
        if (lambda_nm < 510.0)                      return { 30, 170, 170};
        if (lambda_nm < 580.0)                      return { 50, 160,  60};
        if (lambda_nm < 645.0)                      return {220, 170,  30};
        return                                             {200,  60,  60};
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


    // -----------------------------------------------------------------------
    // Thermal colour ramp on light theme.
    //
    //   Starts at a calm light blue (ambient ~ 300 K) and progresses through
    //   green / yellow / orange / red towards bright red (thermal runaway).
    //   Designed to read intuitively on a white-ish backdrop.
    // -----------------------------------------------------------------------
    [[nodiscard]] inline sf::Color thermalColor(float t) noexcept {
        if (t <= 0.0f) return {180, 200, 230};
        if (t >= 1.0f) return {200,  30,  30};

        struct Stop { float t; sf::Color c; };
        static const Stop stops[] = {
            {0.00f, sf::Color(180, 200, 230)},  // light blue (cool)
            {0.20f, sf::Color(120, 200, 220)},  // teal
            {0.40f, sf::Color(120, 200, 100)},  // green
            {0.60f, sf::Color(240, 220,  60)},  // yellow
            {0.80f, sf::Color(240, 130,  40)},  // orange
            {0.95f, sf::Color(220,  60,  40)},  // red
            {1.00f, sf::Color(200,  30,  30)},  // dark red (runaway)
        };
        constexpr int N = sizeof(stops) / sizeof(stops[0]);
        for (int k = 0; k + 1 < N; ++k) {
            if (t <= stops[k + 1].t) {
                const float u = (t - stops[k].t)
                              / (stops[k + 1].t - stops[k].t);
                const auto a = stops[k].c, b = stops[k + 1].c;
                return sf::Color(
                    static_cast<std::uint8_t>(a.r + u * (b.r - a.r)),
                    static_cast<std::uint8_t>(a.g + u * (b.g - a.g)),
                    static_cast<std::uint8_t>(a.b + u * (b.b - a.b)));
            }
        }
        return {200, 30, 30};
    }


    // -----------------------------------------------------------------------
    // Carrier-density heatmap (perceptually-uniform-ish, blue->yellow->red).
    // Replaces the dark viridis ramp so cells stay legible on white.
    // -----------------------------------------------------------------------
    [[nodiscard]] inline sf::Color heatColor(float t) noexcept {
        if (t <= 0.0f) return {220, 230, 245};
        if (t >= 1.0f) return {200,  30,  30};

        struct Stop { float t; sf::Color c; };
        static const Stop stops[] = {
            {0.00f, sf::Color(220, 230, 245)},  // pale blue
            {0.20f, sf::Color(140, 180, 230)},  // sky blue
            {0.40f, sf::Color( 60, 160, 200)},  // teal
            {0.60f, sf::Color( 90, 180,  90)},  // green
            {0.80f, sf::Color(240, 200,  40)},  // gold
            {1.00f, sf::Color(200,  30,  30)},  // red
        };
        constexpr int N = sizeof(stops) / sizeof(stops[0]);
        for (int k = 0; k + 1 < N; ++k) {
            if (t <= stops[k + 1].t) {
                const float u = (t - stops[k].t)
                              / (stops[k + 1].t - stops[k].t);
                const auto a = stops[k].c, b = stops[k + 1].c;
                return sf::Color(
                    static_cast<std::uint8_t>(a.r + u * (b.r - a.r)),
                    static_cast<std::uint8_t>(a.g + u * (b.g - a.g)),
                    static_cast<std::uint8_t>(a.b + u * (b.b - a.b)));
            }
        }
        return {200, 30, 30};
    }

} // namespace palette
