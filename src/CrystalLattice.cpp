#include "CrystalLattice.hpp"

#include <algorithm>
#include <cmath>

// =============================================================================
// Construction / geometry
// =============================================================================
CrystalLattice::CrystalLattice(const sf::Vector2f& topLeft,
                               const sf::Vector2f& size,
                               int rows,
                               int cols)
    : m_topLeft(topLeft), m_size(size)
    , m_rows(rows), m_cols(cols)
    , m_cellW(size.x / static_cast<float>(cols))
    , m_cellH(size.y / static_cast<float>(rows))
    , m_rng(std::random_device{}())
{
    m_atoms.reserve(static_cast<std::size_t>(m_rows * m_cols));
    for (int r = 0; r < m_rows; ++r) {
        for (int c = 0; c < m_cols; ++c) {
            Atom a;
            a.pos = { m_topLeft.x + m_cellW * (c + 0.5f),
                      m_topLeft.y + m_cellH * (r + 0.5f) };
            a.type = AtomType::Silicon;
            m_atoms.push_back(a);
        }
    }
}

// =============================================================================
// Log-scaled carrier-count mapping: physical densities span many orders of
// magnitude (1e10..1e20), so we squash them into a bounded visual count.
// =============================================================================
int CrystalLattice::targetCarrierCount(double concentration) const {
    if (concentration <= 0.0) return 0;
    const double logN = std::log10(concentration);
    const double scaled = (logN - 10.0) * 8.0;
    return static_cast<int>(std::clamp(scaled, 0.0, 120.0));
}

// =============================================================================
// Dopant placement + carrier resample.
// =============================================================================
void CrystalLattice::rebuild(const PhysicsEngine& physics) {
    // Reset to pure Si.
    for (auto& a : m_atoms) a.type = AtomType::Silicon;

    const DopingType dt = physics.getDopingType();

    if (dt != DopingType::Intrinsic && physics.getDopingConcentration() > 0.0) {
        const double logN = std::log10(physics.getDopingConcentration());
        const double frac = std::clamp((logN - 14.0) / 6.0 * 0.23 + 0.02,
                                       0.0, 0.30);
        const int total  = static_cast<int>(m_atoms.size());
        const int numDop = static_cast<int>(total * frac);

        std::uniform_int_distribution<int> pick(0, total - 1);
        const AtomType target = (dt == DopingType::NType)
                              ? AtomType::Phosphorus
                              : AtomType::Boron;
        for (int i = 0; i < numDop; ++i) {
            m_atoms[pick(m_rng)].type = target;
        }
    }

    resampleCarriers(physics);
}

// =============================================================================
// Populate the carrier cloud.
//
// Thermal carriers reflect n and p from the physics engine.  Optical
// carriers (marked with `optical = true`) reflect Delta n from photon
// absorption and are drawn in a brighter / cyan-ish colour so the user
// can distinguish "dark-equilibrium" carriers from "photo-generated" ones.
// =============================================================================
void CrystalLattice::resampleCarriers(const PhysicsEngine& physics) {
    m_carriers.clear();

    const int nElectrons     = targetCarrierCount(physics.getElectronConcentration());
    const int nHoles         = targetCarrierCount(physics.getHoleConcentration());
    const int nOpticalPairs  = targetCarrierCount(physics.getExcessCarrierDensity());

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

    // Thermal carriers
    spawn(nElectrons, true,  false);
    spawn(nHoles,     false, false);
    // Optical electron-hole pairs (generated in equal number)
    spawn(nOpticalPairs, true,  true);
    spawn(nOpticalPairs, false, true);
}

// =============================================================================
// Per-frame update.
//
//   * Brownian kick  -> random-walk appearance.
//   * Velocity damping keeps the mean speed bounded.
//   * Lorentz force  F = q (v x B)  with B perpendicular to the screen
//     manifests as a per-frame rotation of every carrier's velocity vector.
//     For a 2D (x,y) velocity and B along +z, v x B = (v_y B, -v_x B, 0),
//     so the Newton equation m dv/dt = q(v x B) is simply the cyclotron
//     motion with angular frequency  omega_c = q B / m*  (see Kasap 2.4,
//     Kittel Ch. 8).  Electrons (q<0) and holes (q>0) therefore orbit in
//     opposite senses -- which is exactly the visual signature of the Hall
//     effect.  We pick a visualization-friendly omega scale rather than the
//     literal m*, since the simulator is pedagogical, not quantitative, in
//     space-time coordinates.
//   * Wall bounces preserve carriers inside the visualization rectangle.
// =============================================================================
void CrystalLattice::update(float dt, const PhysicsEngine& physics) {
    std::uniform_real_distribution<float> kick(-40.f, 40.f);

    const float xMin = m_topLeft.x + 2.f;
    const float xMax = m_topLeft.x + m_size.x - 2.f;
    const float yMin = m_topLeft.y + 2.f;
    const float yMax = m_topLeft.y + m_size.y - 2.f;

    // Visualization-scale cyclotron frequency [rad/s per Tesla].
    constexpr float kLorentzScale = 3.0f;
    const float B     = static_cast<float>(physics.getMagneticField());
    const float omega = kLorentzScale * B;   // |omega_c| scaled for screen
    const bool  magnetic = (B != 0.0f);

    // Cache rotation matrix components if omega * dt is constant this frame.
    const float dtheta = omega * dt;
    const float cs     = magnetic ? std::cos(dtheta) : 1.f;
    const float sn     = magnetic ? std::sin(dtheta) : 0.f;

    for (auto& c : m_carriers) {
        // --- Thermal random walk --------------------------------------------
        c.vel.x += kick(m_rng) * dt;
        c.vel.y += kick(m_rng) * dt;
        c.vel   *= 0.995f;

        // --- Lorentz rotation (opposite sense for electrons vs holes) -------
        if (magnetic) {
            const float sgn = c.electron ? +1.f : -1.f;
            const float vx_new = c.vel.x * cs - sgn * c.vel.y * sn;
            const float vy_new = sgn * c.vel.x * sn + c.vel.y * cs;
            c.vel = { vx_new, vy_new };
        }

        c.pos += c.vel * dt;

        // --- Wall bouncing --------------------------------------------------
        if (c.pos.x < xMin) { c.pos.x = xMin; c.vel.x = -c.vel.x; }
        if (c.pos.x > xMax) { c.pos.x = xMax; c.vel.x = -c.vel.x; }
        if (c.pos.y < yMin) { c.pos.y = yMin; c.vel.y = -c.vel.y; }
        if (c.pos.y > yMax) { c.pos.y = yMax; c.vel.y = -c.vel.y; }
    }
}

// =============================================================================
// Rendering
// =============================================================================
void CrystalLattice::draw(sf::RenderTarget& target) const {
    // Background
    sf::RectangleShape bg(m_size);
    bg.setPosition(m_topLeft);
    bg.setFillColor(sf::Color(18, 22, 30));
    bg.setOutlineColor(sf::Color(80, 90, 110));
    bg.setOutlineThickness(1.f);
    target.draw(bg);

    // Bonds
    sf::VertexArray bonds(sf::Lines);
    for (int r = 0; r < m_rows; ++r) {
        for (int c = 0; c < m_cols; ++c) {
            const auto& a = m_atoms[r * m_cols + c];
            const sf::Color col(70, 80, 100);
            if (c + 1 < m_cols) {
                const auto& b = m_atoms[r * m_cols + (c + 1)];
                bonds.append(sf::Vertex(a.pos, col));
                bonds.append(sf::Vertex(b.pos, col));
            }
            if (r + 1 < m_rows) {
                const auto& b = m_atoms[(r + 1) * m_cols + c];
                bonds.append(sf::Vertex(a.pos, col));
                bonds.append(sf::Vertex(b.pos, col));
            }
        }
    }
    target.draw(bonds);

    // Atoms
    const float radius = std::min(m_cellW, m_cellH) * 0.25f;
    sf::CircleShape atomShape(radius);
    atomShape.setOrigin(radius, radius);
    for (const auto& a : m_atoms) {
        switch (a.type) {
            case AtomType::Silicon:
                atomShape.setFillColor(sf::Color(120, 170, 220));
                atomShape.setOutlineColor(sf::Color(200, 220, 240));
                atomShape.setOutlineThickness(1.f);
                break;
            case AtomType::Phosphorus:
                atomShape.setFillColor(sf::Color( 80, 220, 120));
                atomShape.setOutlineColor(sf::Color(220, 255, 200));
                atomShape.setOutlineThickness(1.5f);
                break;
            case AtomType::Boron:
                atomShape.setFillColor(sf::Color(230, 110, 110));
                atomShape.setOutlineColor(sf::Color(255, 200, 200));
                atomShape.setOutlineThickness(1.5f);
                break;
        }
        atomShape.setPosition(a.pos);
        target.draw(atomShape);
    }

    // Carriers.  Optical carriers are rendered slightly larger and with a
    // halo so illumination is visually obvious.
    sf::CircleShape dot(3.f);
    dot.setOrigin(3.f, 3.f);
    for (const auto& c : m_carriers) {
        sf::Color fill, outline;
        if (c.electron) {
            fill    = c.optical ? sf::Color(140, 255, 255)
                                : sf::Color(255, 230,  80);
            outline = sf::Color(255, 255, 220);
        } else {
            fill    = c.optical ? sf::Color(255, 180, 255)
                                : sf::Color(255,  80, 200);
            outline = sf::Color(255, 220, 240);
        }
        dot.setFillColor(fill);
        dot.setOutlineColor(outline);
        dot.setOutlineThickness(c.optical ? 1.5f : 1.f);
        dot.setRadius(c.optical ? 4.f : 3.f);
        dot.setOrigin(c.optical ? 4.f : 3.f, c.optical ? 4.f : 3.f);
        dot.setPosition(c.pos);
        target.draw(dot);
    }
}
