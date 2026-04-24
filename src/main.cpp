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

#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include "BandDiagram.hpp"
#include "CrystalLattice.hpp"
#include "PhysicsEngine.hpp"
#include "UIPanel.hpp"


// -----------------------------------------------------------------------------
// Font loader: tries a few common locations so the project works out-of-the-box
// on Linux and Windows. If no font can be loaded, the program still runs but
// text labels will be missing.
// -----------------------------------------------------------------------------
static bool loadFont(sf::Font& font) {
    const std::vector<std::string> candidates = {
        // Project-local fallback (see assets/README.md).
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

    for (const auto& path : candidates) {
        if (font.loadFromFile(path)) {
            std::cout << "[info] Loaded font: " << path << "\n";
            return true;
        }
    }
    std::cerr << "[warn] Could not load any font. "
                 "Place a TTF at assets/font.ttf.\n";
    return false;
}


int main() {
    // -------------------------------------------------------------------------
    // Window
    // -------------------------------------------------------------------------
    constexpr unsigned winW = 1400;
    constexpr unsigned winH = 900;

    sf::RenderWindow window(
        sf::VideoMode(winW, winH),
        "Interactive Semiconductor Bandgap & Doping Simulator",
        sf::Style::Titlebar | sf::Style::Close);
    window.setFramerateLimit(60);

    // -------------------------------------------------------------------------
    // Resources
    // -------------------------------------------------------------------------
    sf::Font font;
    loadFont(font);   // continue even if this fails

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

        // ---- Event handling ------------------------------------------------
        sf::Event event;
        while (window.pollEvent(event)) {
            if (event.type == sf::Event::Closed)
                window.close();
            if (event.type == sf::Event::KeyPressed &&
                event.key.code == sf::Keyboard::Escape)
                window.close();

            if (ui->handleEvent(event, mousePos, *physics)) {
                lattice->rebuild(*physics);
            }
        }

        // ---- Update --------------------------------------------------------
        const float dt = clock.restart().asSeconds();
        lattice->update(dt, *physics);

        // ---- Render --------------------------------------------------------
        window.clear(sf::Color(10, 12, 18));

        lattice->draw(window);
        bands->draw(window, *physics, font);
        ui->draw(window, *physics, font);

        window.display();
    }

    return 0;
}
