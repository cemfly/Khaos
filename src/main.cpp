// =============================================================================
// Interactive Semiconductor Bandgap & Doping Simulator
// -----------------------------------------------------------------------------
// Entry point. Creates the SFML window, owns the physics engine and the three
// visualization components (crystal lattice, band diagram, UI panel), and
// drives the event / update / render loop.
//
// Layout (1400 x 900 window):
//
//   +-----------------------------+-----------------------------+
//   |                             |                             |
//   |   Crystal lattice (2D)      |   Energy band diagram       |
//   |   - Si atoms                |   + Fermi-Dirac curve       |
//   |   - Dopant atoms (P / B)    |                             |
//   |   - Free electrons / holes  |                             |
//   |                             |                             |
//   +-----------------------------+-----------------------------+
//   |   Control panel: Temperature, Doping concentration,       |
//   |                  Intrinsic / n-type / p-type selection    |
//   +-----------------------------------------------------------+
// =============================================================================

#include <SFML/Graphics.hpp>

#include <array>
#include <filesystem>
#include <iostream>
#include <memory>
#include <optional>
#include <string_view>

#include "BandDiagram.hpp"
#include "CrystalLattice.hpp"
#include "Palette.hpp"
#include "PhysicsEngine.hpp"
#include "UIPanel.hpp"


// -----------------------------------------------------------------------------
// Font loader
// -----------------------------------------------------------------------------
// Probes a fixed list of well-known font locations (cross-platform). Only the
// paths that actually exist on disk are reported to sf::Font::openFromFile so
// that the console is not polluted with SFML's "file not found" warnings for
// paths that are obviously not going to exist (e.g. Linux paths on Windows).
// Returns true on the first successful load.
// -----------------------------------------------------------------------------
namespace {

constexpr std::array<std::string_view, 7> kFontCandidates{
    // Project-local (preferred -- copied into build/ by CMake).
    "assets/font.ttf",
    "../assets/font.ttf",
    // Common Linux locations.
    "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf",
    "/usr/share/fonts/truetype/liberation/LiberationSans-Regular.ttf",
    "/usr/share/fonts/TTF/DejaVuSans.ttf",
    // Common Windows locations.
    "C:/Windows/Fonts/arial.ttf",
    "C:/Windows/Fonts/segoeui.ttf",
};

[[nodiscard]] bool loadFont(sf::Font& font) {
    for (const auto path : kFontCandidates) {
        const std::filesystem::path fp{path};
        std::error_code ec;
        if (!std::filesystem::exists(fp, ec)) continue;
        if (font.openFromFile(fp)) {
            std::cout << "[info] Loaded font: " << fp.string() << '\n';
            return true;
        }
    }
    std::cerr << "[warn] No usable font found. "
                 "Drop a .ttf at assets/font.ttf.\n";
    return false;
}

} // namespace


int main() {
    // -------------------------------------------------------------------------
    // Window
    // -------------------------------------------------------------------------
    constexpr unsigned winW = 1400;
    constexpr unsigned winH = 900;

    sf::RenderWindow window(
        sf::VideoMode({winW, winH}),
        "Interactive Semiconductor Bandgap & Doping Simulator",
        sf::Style::Titlebar | sf::Style::Close);
    window.setFramerateLimit(60);

    // -------------------------------------------------------------------------
    // Resources
    // -------------------------------------------------------------------------
    sf::Font font;
    (void)loadFont(font);   // non-fatal if it fails

    // -------------------------------------------------------------------------
    // Simulation components  (owned by smart pointers; created once here).
    // -------------------------------------------------------------------------
    auto physics  = std::make_unique<PhysicsEngine>();

    // Phase 3 layout.  The UI panel is taller (260 px) to host the
    // wavelength + B-field sliders and the extra mode-toggle buttons, so
    // the lattice and band-diagram panels are shrunk vertically to fit.
    auto lattice  = std::make_unique<CrystalLattice>(
        sf::Vector2f{20.f, 20.f},          // top-left
        sf::Vector2f{640.f, 580.f},        // size
        18,                                 // rows
        18);                                // cols

    auto bands    = std::make_unique<BandDiagram>(
        sf::Vector2f{680.f, 20.f},
        sf::Vector2f{700.f, 580.f});

    auto ui       = std::make_unique<UIPanel>(
        sf::Vector2f{20.f, 620.f},
        sf::Vector2f{1360.f, 260.f});

    // Initial build of the lattice visualization.
    lattice->rebuild(*physics);

    // -------------------------------------------------------------------------
    // Main loop
    // -------------------------------------------------------------------------
    sf::Clock clock;
    while (window.isOpen()) {
        const sf::Vector2f mousePos =
            window.mapPixelToCoords(sf::Mouse::getPosition(window));

        // ---- Event handling (SFML 3 optional-based polling) ---------------
        while (const std::optional<sf::Event> event = window.pollEvent()) {
            if (event->is<sf::Event::Closed>()) {
                window.close();
            }
            else if (const auto* keyEv =
                         event->getIf<sf::Event::KeyPressed>())
            {
                if (keyEv->code == sf::Keyboard::Key::Escape)
                    window.close();
            }

            if (ui->handleEvent(*event, mousePos, *physics)) {
                lattice->rebuild(*physics);
            }
        }

        // ---- Update --------------------------------------------------------
        const float dt = clock.restart().asSeconds();
        lattice->update(dt, *physics);

        // ---- Render --------------------------------------------------------
        window.clear(palette::WindowBg);

        lattice->draw(window);
        bands->draw(window, *physics, font);
        ui->draw(window, *physics, font);

        window.display();
    }

    return 0;
}
