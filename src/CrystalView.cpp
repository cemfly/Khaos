// =============================================================================
// CrystalView.cpp
//
//   Author : dex / cemfly-april2026
//   License: MIT
// =============================================================================

#include "CrystalView.hpp"

#include <algorithm>
#include <cmath>

#include "Palette.hpp"


// =============================================================================
// Construction
// =============================================================================
CrystalView::CrystalView(unsigned int textureSize)
    : m_topLeft{0.f, 0.f}
    , m_size{static_cast<float>(textureSize), static_cast<float>(textureSize)}
    , m_rows(20)
    , m_cols(20)
    , m_cellW(m_size.x / static_cast<float>(m_cols))
    , m_cellH(m_size.y / static_cast<float>(m_rows))
    , m_rng(std::random_device{}())
{
    if (!m_rt.resize({textureSize, textureSize})) {
        // Texture creation failed; renderTexture().getSize() will be {0,0}.
    }
    m_rt.setSmooth(true);

    m_atoms.reserve(static_cast<std::size_t>(m_rows * m_cols));
    for (int r = 0; r < m_rows; ++r) {
        for (int c = 0; c < m_cols; ++c) {
            Atom a;
            a.pos = {
                m_topLeft.x + m_cellW * (c + 0.5f),
                m_topLeft.y + m_cellH * (r + 0.5f)
            };
            a.type = AtomType::Host;
            m_atoms.push_back(a);
        }
    }
}


// =============================================================================
// Carrier-count mapping (log-scaled visual density)
// =============================================================================
int CrystalView::targetCarrierCount(double concentration) const noexcept {
    if (concentration <= 0.0) return 0;
    const double logN   = std::log10(concentration);
    const double scaled = (logN - 8.0) * 6.0;
    return static_cast<int>(std::clamp(scaled, 0.0, 120.0));
}


// =============================================================================
// Rebuild atoms + carriers
// =============================================================================
void CrystalView::rebuild(const PhysicsEngine& physics) {
    for (auto& a : m_atoms) a.type = AtomType::Host;

    const auto dt = physics.getDopingType();
    if (dt != DopingType::Intrinsic && physics.getDopingConcentration() > 0.0) {
        const double logN = std::log10(physics.getDopingConcentration());
        const double frac = std::clamp((logN - 14.0) / 6.0 * 0.23 + 0.02,
                                       0.0, 0.30);
        const int total  = static_cast<int>(m_atoms.size());
        const int numDop = static_cast<int>(total * frac);

        std::uniform_int_distribution<int> pick(0, total - 1);
        const AtomType target = (dt == DopingType::NType) ? AtomType::Donor
                                                          : AtomType::Acceptor;
        for (int i = 0; i < numDop; ++i)
            m_atoms[pick(m_rng)].type = target;
    }

    resampleCarriers(physics);
}


void CrystalView::resampleCarriers(const PhysicsEngine& physics) {
    m_carriers.clear();

    const int nElectrons    = targetCarrierCount(physics.getElectronConcentration());
    const int nHoles        = targetCarrierCount(physics.getHoleConcentration());
    const int nOpticalPairs = targetCarrierCount(physics.getExcessCarrierDensity());

    std::uniform_real_distribution<float> xDist(m_topLeft.x + 4.f,
                                                m_topLeft.x + m_size.x - 4.f);
    std::uniform_real_distribution<float> yDist(m_topLeft.y + 4.f,
                                                m_topLeft.y + m_size.y - 4.f);
    std::uniform_real_distribution<float> vDist(-25.f, 25.f);

    auto spawn = [&](int count, bool electron, bool optical) {
        for (int i = 0; i < count; ++i) {
            Carrier c;
            c.pos      = { xDist(m_rng), yDist(m_rng) };
            c.vel      = { vDist(m_rng), vDist(m_rng) };
            c.electron = electron;
            c.optical  = optical;
            m_carriers.push_back(c);
        }
    };

    spawn(nElectrons, true,  false);
    spawn(nHoles,     false, false);
    spawn(nOpticalPairs, true,  true);
    spawn(nOpticalPairs, false, true);
}


// =============================================================================
// Per-frame carrier update
// =============================================================================
void CrystalView::update(float dt, const PhysicsEngine& physics) {
    std::uniform_real_distribution<float> kick(-40.f, 40.f);

    const float xMin = m_topLeft.x + 2.f;
    const float xMax = m_topLeft.x + m_size.x - 2.f;
    const float yMin = m_topLeft.y + 2.f;
    const float yMax = m_topLeft.y + m_size.y - 2.f;

    constexpr float kLorentzScale = 3.0f;
    const float B        = static_cast<float>(physics.getMagneticField());
    const float omega    = kLorentzScale * B;
    const bool  magnetic = (B != 0.0f);

    const float dtheta = omega * dt;
    const float cs     = magnetic ? std::cos(dtheta) : 1.f;
    const float sn     = magnetic ? std::sin(dtheta) : 0.f;

    for (auto& c : m_carriers) {
        c.vel.x += kick(m_rng) * dt;
        c.vel.y += kick(m_rng) * dt;
        c.vel   *= 0.995f;

        if (magnetic) {
            const float sgn   = c.electron ? +1.f : -1.f;
            const float vx_n  = c.vel.x * cs - sgn * c.vel.y * sn;
            const float vy_n  = sgn * c.vel.x * sn + c.vel.y * cs;
            c.vel = { vx_n, vy_n };
        }

        c.pos += c.vel * dt;

        if (c.pos.x < xMin) { c.pos.x = xMin; c.vel.x = -c.vel.x; }
        if (c.pos.x > xMax) { c.pos.x = xMax; c.vel.x = -c.vel.x; }
        if (c.pos.y < yMin) { c.pos.y = yMin; c.vel.y = -c.vel.y; }
        if (c.pos.y > yMax) { c.pos.y = yMax; c.vel.y = -c.vel.y; }
    }
}


// =============================================================================
// Carrier heatmap (n(x,y) from DriftDiffusion grid)
// =============================================================================
void CrystalView::drawCarrierHeatmap(const DriftDiffusion& dd) {
    const float maxv = std::max(dd.maxValue(), 1.0e-3f);
    const int   W    = dd.width();
    const int   H    = dd.height();
    const float cw   = m_size.x / static_cast<float>(W);
    const float ch   = m_size.y / static_cast<float>(H);

    sf::VertexArray quads(sf::PrimitiveType::Triangles);
    quads.resize(static_cast<std::size_t>(W * H * 6));

    std::size_t v = 0;
    for (int j = 0; j < H; ++j) {
        for (int i = 0; i < W; ++i) {
            const float n  = dd.at(i, j);
            const float t  = std::clamp(n / maxv, 0.0f, 1.0f);
            sf::Color c    = palette::heatColor(t);
            c.a            = static_cast<std::uint8_t>(160.0f * t);

            const float x0 = m_topLeft.x + i * cw;
            const float y0 = m_topLeft.y + j * ch;
            const float x1 = x0 + cw;
            const float y1 = y0 + ch;

            quads[v++] = sf::Vertex{{x0, y0}, c};
            quads[v++] = sf::Vertex{{x1, y0}, c};
            quads[v++] = sf::Vertex{{x1, y1}, c};
            quads[v++] = sf::Vertex{{x0, y0}, c};
            quads[v++] = sf::Vertex{{x1, y1}, c};
            quads[v++] = sf::Vertex{{x0, y1}, c};
        }
    }
    m_rt.draw(quads);
}


// =============================================================================
// Thermal heatmap (T(x,y))
//
// Mapping:  300 K (ambient)  ->  blue
//           600 K            ->  yellow / orange
//           >900 K           ->  red / white-hot (thermal runaway)
// =============================================================================
void CrystalView::drawThermalHeatmap(const DriftDiffusion& dd) {
    const float Tmin = 300.0f;
    const float Tmax = std::max(dd.maxTemperature(), Tmin + 1.0f);
    const float Tspan = std::max(Tmax - Tmin, 1.0f);

    const int   W = dd.width();
    const int   H = dd.height();
    const float cw = m_size.x / static_cast<float>(W);
    const float ch = m_size.y / static_cast<float>(H);

    sf::VertexArray quads(sf::PrimitiveType::Triangles);
    quads.resize(static_cast<std::size_t>(W * H * 6));

    std::size_t v = 0;
    for (int j = 0; j < H; ++j) {
        for (int i = 0; i < W; ++i) {
            const float T  = dd.temperatureAt(i, j);
            const float t  = std::clamp((T - Tmin) / Tspan, 0.0f, 1.0f);
            sf::Color c    = palette::thermalColor(t);
            // Stronger alpha than carrier heatmap so hot regions glow
            // visibly through the lattice.
            c.a            = static_cast<std::uint8_t>(60.0f + 140.0f * t);

            const float x0 = m_topLeft.x + i * cw;
            const float y0 = m_topLeft.y + j * ch;
            const float x1 = x0 + cw;
            const float y1 = y0 + ch;

            quads[v++] = sf::Vertex{{x0, y0}, c};
            quads[v++] = sf::Vertex{{x1, y0}, c};
            quads[v++] = sf::Vertex{{x1, y1}, c};
            quads[v++] = sf::Vertex{{x0, y0}, c};
            quads[v++] = sf::Vertex{{x1, y1}, c};
            quads[v++] = sf::Vertex{{x0, y1}, c};
        }
    }
    m_rt.draw(quads);
}


// =============================================================================
// BJT region overlay (Emitter / Base / Collector)
// =============================================================================
void CrystalView::drawBjtRegions(const DriftDiffusion& dd) {
    if (dd.deviceMode() != DeviceMode::NpnBjt) return;

    const int   W  = dd.width();
    const int   H  = dd.height();
    const float cw = m_size.x / static_cast<float>(W);
    const float ch = m_size.y / static_cast<float>(H);

    sf::VertexArray quads(sf::PrimitiveType::Triangles);
    quads.resize(static_cast<std::size_t>(W * H * 6));

    std::size_t v = 0;
    for (int j = 0; j < H; ++j) {
        for (int i = 0; i < W; ++i) {
            const auto r = dd.regionAt(i, j);
            sf::Color c{0, 0, 0, 0};
            switch (r) {
                case CellRegion::Emitter:
                    c = sf::Color(80, 180, 255,  60);  // blue tint
                    break;
                case CellRegion::Base:
                    c = sf::Color(255, 130, 130,  90); // pink tint -- the
                                                       // narrow active region
                    break;
                case CellRegion::Collector:
                    c = sf::Color(140, 255, 180,  50); // green tint
                    break;
                case CellRegion::Bulk:
                default:
                    continue;
            }

            const float x0 = m_topLeft.x + i * cw;
            const float y0 = m_topLeft.y + j * ch;
            const float x1 = x0 + cw;
            const float y1 = y0 + ch;

            quads[v++] = sf::Vertex{{x0, y0}, c};
            quads[v++] = sf::Vertex{{x1, y0}, c};
            quads[v++] = sf::Vertex{{x1, y1}, c};
            quads[v++] = sf::Vertex{{x0, y0}, c};
            quads[v++] = sf::Vertex{{x1, y1}, c};
            quads[v++] = sf::Vertex{{x0, y1}, c};
        }
    }
    quads.resize(v);
    m_rt.draw(quads);
}


// =============================================================================
// Lattice (atoms + bonds)
// =============================================================================
void CrystalView::drawLattice(const PhysicsEngine& physics) {
    sf::VertexArray bonds(sf::PrimitiveType::Lines);
    for (int r = 0; r < m_rows; ++r) {
        for (int c = 0; c < m_cols; ++c) {
            const auto& a = m_atoms[r * m_cols + c];
            const sf::Color col(160, 175, 200);   // light blue-grey on white
            if (c + 1 < m_cols) {
                const auto& b = m_atoms[r * m_cols + (c + 1)];
                bonds.append(sf::Vertex{a.pos, col});
                bonds.append(sf::Vertex{b.pos, col});
            }
            if (r + 1 < m_rows) {
                const auto& b = m_atoms[(r + 1) * m_cols + c];
                bonds.append(sf::Vertex{a.pos, col});
                bonds.append(sf::Vertex{b.pos, col});
            }
        }
    }
    m_rt.draw(bonds);

    const auto& mat = physics.getMaterial();
    const sf::Color hostCol(mat.atomR, mat.atomG, mat.atomB);

    const float radius = std::min(m_cellW, m_cellH) * 0.25f;
    sf::CircleShape atom(radius);
    atom.setOrigin({radius, radius});
    for (const auto& a : m_atoms) {
        switch (a.type) {
            case AtomType::Host:
                atom.setFillColor(hostCol);
                atom.setOutlineColor(sf::Color(40, 50, 70, 140));
                atom.setOutlineThickness(1.f);
                break;
            case AtomType::Donor:
                atom.setFillColor(palette::Phosphorus);
                atom.setOutlineColor(sf::Color(20, 70, 30, 180));
                atom.setOutlineThickness(1.5f);
                break;
            case AtomType::Acceptor:
                atom.setFillColor(palette::Boron);
                atom.setOutlineColor(sf::Color(100, 30, 30, 180));
                atom.setOutlineThickness(1.5f);
                break;
        }
        atom.setPosition(a.pos);
        m_rt.draw(atom);
    }
}


// =============================================================================
// Carriers
// =============================================================================
void CrystalView::drawCarriers() {
    sf::CircleShape dot(3.f);
    dot.setOrigin({3.f, 3.f});
    for (const auto& c : m_carriers) {
        const sf::Color fill = c.electron
            ? (c.optical ? palette::ElectronOpt : palette::Electron)
            : (c.optical ? palette::HoleOpt     : palette::Hole);
        // Dark thin outlines so dots stay legible on the light background.
        const sf::Color outline = c.electron ? sf::Color(110, 60, 20, 200)
                                             : sf::Color( 60, 20, 80, 200);
        dot.setFillColor(fill);
        dot.setOutlineColor(outline);
        dot.setOutlineThickness(c.optical ? 1.5f : 1.f);
        const float r = c.optical ? 4.f : 3.f;
        dot.setRadius(r);
        dot.setOrigin({r, r});
        dot.setPosition(c.pos);
        m_rt.draw(dot);
    }
}


// =============================================================================
// Lorentz vector field
// =============================================================================
void CrystalView::drawVectorField(const PhysicsEngine& physics) {
    const float B = static_cast<float>(physics.getMagneticField());
    if (std::abs(B) < 1.0e-6f) return;

    constexpr int   step       = 60;
    const float     maxLen     = step * 0.45f;
    const float     normalised = std::min(std::abs(B) / 10.0f, 1.0f);
    const float     len        = maxLen * normalised;
    const float     sign       = (B > 0.0f) ? +1.f : -1.f;

    sf::VertexArray lines(sf::PrimitiveType::Lines);
    sf::Color col(palette::Hall.r, palette::Hall.g, palette::Hall.b, 200);

    for (float y = step / 2.f; y < m_size.y; y += step) {
        for (float x = step / 2.f; x < m_size.x; x += step) {
            const float x0 = m_topLeft.x + x;
            const float y0 = m_topLeft.y + y;
            const float yEnd = y0 + sign * len;

            lines.append(sf::Vertex{{x0, y0},   col});
            lines.append(sf::Vertex{{x0, yEnd}, col});

            lines.append(sf::Vertex{{x0,        yEnd},                 col});
            lines.append(sf::Vertex{{x0 - 4.f, yEnd - sign * 6.f},     col});
            lines.append(sf::Vertex{{x0,        yEnd},                 col});
            lines.append(sf::Vertex{{x0 + 4.f, yEnd - sign * 6.f},     col});
        }
    }
    m_rt.draw(lines);
}


// =============================================================================
// Master compositor
// =============================================================================
void CrystalView::render(const PhysicsEngine&  physics,
                         const DriftDiffusion& dd,
                         HeatmapMode           heatmapMode,
                         bool                  showVectorField)
{
    m_rt.clear(palette::ViewBg);

    switch (heatmapMode) {
        case HeatmapMode::Carriers: drawCarrierHeatmap(dd); break;
        case HeatmapMode::Thermal:  drawThermalHeatmap(dd); break;
        case HeatmapMode::None:
        default: break;
    }

    drawBjtRegions(dd);
    drawLattice(physics);
    drawCarriers();
    if (showVectorField) drawVectorField(physics);

    m_rt.display();
}
