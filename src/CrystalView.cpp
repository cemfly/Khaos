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
// Per-frame carrier transport -- coupled to the real DriftDiffusion field
// -----------------------------------------------------------------------------
// Three rules drive the visualisation, designed so the eye reads "current"
// rather than "bouncing balls":
//
//   1. Absorbing + injecting X-boundaries (Continuous Flow)
//      A carrier that crosses the left or right wall is *not* reflected;
//      it is teleported to the opposite side at a random y. This emulates
//      the steady-state current loop: a particle leaving the right N-contact
//      is replaced by an equivalent particle entering from the left
//      P-contact (charge neutrality of the external circuit).
//
//   2. Depletion-region culling (Junction empties out)
//      A visual electron in a cell where the *real* DD density n(x,y) is
//      below ~3 n_i is in the depletion region; it gets relocated to a
//      randomly-sampled N-majority bulk cell. Holes likewise relocate to a
//      P-majority cell. The junction stays visually empty, which is what
//      the physics says (Sze Sec. 2.2: free carriers are swept out of the
//      depletion region).
//
//   3. Drift > diffusion (Visual scaling)
//      The Einstein thermal sigma is damped by VISUAL_DIFFUSION_DAMPER
//      (= 0.05) so the deterministic drift dominates the rendered motion.
//      Without it the Brownian kicks would mask the v_sat sweep-out.
//
// References: Sze Sec. 1.4 (drift), 2.2 (depletion), Pierret Sec. 3.4
// (Einstein), Caughey-Thomas IEEE Proc. 55 (1967) for v_sat.
// Zero-allocation: m_carriers iterated in-place; all distributions are
// stack-local stateless objects; m_rng reused.
// =============================================================================
void CrystalView::update(float dt,
                         const PhysicsEngine&  physics,
                         const DriftDiffusion& dd)
{
    if (dt <= 0.0f || m_carriers.empty()) return;

    // ---- Grid / pixel scaling -----------------------------------------------
    const int   W = dd.width();
    const int   H = dd.height();
    if (W <= 0 || H <= 0) return;
    const float pitch_cm  = std::max(dd.cellPitchCm(), 1.0e-9f);
    const float px_per_cm = m_size.x / (static_cast<float>(W) * pitch_cm);
    const float cellW_px  = m_size.x / static_cast<float>(W);
    const float cellH_px  = m_size.y / static_cast<float>(H);

    // ---- Visual-time compression (v_sat -> ~250 px/s on screen) -------------
    const auto& mat = physics.getMaterial();
    const float v_sat_avg = static_cast<float>(
        0.5 * (mat.v_sat_n + mat.v_sat_p));
    constexpr float V_SAT_PIXEL = 250.0f;
    const float tau_vis = (v_sat_avg * px_per_cm > 0.0f)
        ? V_SAT_PIXEL / (v_sat_avg * px_per_cm)
        : 1.0e-10f;
    const float dt_phys = dt * tau_vis;

    // Rule 3: diffusion damper so drift visually wins.
    constexpr float VISUAL_DIFFUSION_DAMPER = 0.05f;

    // ---- Physical state -----------------------------------------------------
    const double T    = physics.getTemperature();
    const double V_T  = PhysicsEngine::thermalVoltage(T);
    const double n_i  = physics.getIntrinsicCarrier();

    // Depletion threshold: cells where n (for electrons) or p (for holes)
    // falls below this are visually emptied. 3 n_i is a robust marker for
    // the depleted core of the junction; deep N or P bulk is many orders
    // higher, so false positives are negligible.
    const double DEPLETION_THRESHOLD = 3.0 * n_i;

    // ---- Lorentz / Hall (rotates drift vector) -----------------------------
    constexpr float kLorentzScale = 3.0f;
    const float B        = static_cast<float>(physics.getMagneticField());
    const bool  magnetic = (B != 0.0f);
    const float dtheta   = kLorentzScale * B * dt;
    const float cs       = magnetic ? std::cos(dtheta) : 1.f;
    const float sn       = magnetic ? std::sin(dtheta) : 0.f;

    // ---- Window in pixel space ---------------------------------------------
    const float xMin = m_topLeft.x + 2.f;
    const float xMax = m_topLeft.x + m_size.x - 2.f;
    const float yMin = m_topLeft.y + 2.f;
    const float yMax = m_topLeft.y + m_size.y - 2.f;
    const float xSpan = xMax - xMin;

    // ---- Stack-local RNG state (no allocation, stateless distributions) ----
    std::normal_distribution<float>  thermal(0.0f, 1.0f);
    std::uniform_real_distribution<float> uniY(yMin, yMax);
    std::uniform_real_distribution<float> uniXedgeOffset(2.f, 0.05f * xSpan);
    std::uniform_real_distribution<float> jitterCell(-0.45f, 0.45f);
    std::uniform_int_distribution<int>    randCellX(0, W - 1);
    std::uniform_int_distribution<int>    randCellY(0, H - 1);
    std::uniform_real_distribution<float> tinyKick(-12.f, 12.f);
    constexpr float E_FIELD_THRESHOLD = 1.0f;   // V/cm

    const bool painter = (dd.deviceMode() == DeviceMode::Painter);

    // ---- Helper: relocate a culled carrier to a random majority-dopant cell.
    //
    // Up to 6 sample attempts (cheap; cells are 60*40 = 2400 typically).
    // If nothing matches, return false so the caller can fall back to a
    // periodic injection at the screen edge (no infinite loop hazard).
    auto teleportToBulk = [&](Carrier& cr) -> bool {
        if (!painter) return false;
        for (int attempt = 0; attempt < 6; ++attempt) {
            const int i = randCellX(m_rng);
            const int j = randCellY(m_rng);
            const double Nd_c = dd.donorAt(i, j);
            const double Na_c = dd.acceptorAt(i, j);
            const bool ok = cr.electron
                ? (Nd_c > Na_c && Nd_c > 1.0e14)
                : (Na_c > Nd_c && Na_c > 1.0e14);
            if (!ok) continue;
            cr.pos.x = m_topLeft.x + (i + 0.5f + jitterCell(m_rng)) * cellW_px;
            cr.pos.y = m_topLeft.y + (j + 0.5f + jitterCell(m_rng)) * cellH_px;
            cr.vel   = {0.0f, 0.0f};
            return true;
        }
        return false;
    };

    // ---- Helper: respawn at the opposite edge (Rule 1 continuous current).
    //
    // Direction of exit determines which edge to re-inject from. The
    // particle's velocity is reset; the next frame's drift evaluator picks
    // up the local field at the new position and accelerates the carrier
    // back into the device.
    auto injectFromOppositeEdge = [&](Carrier& cr, bool exitedRight) {
        cr.pos.x = exitedRight
            ? xMin + uniXedgeOffset(m_rng)
            : xMax - uniXedgeOffset(m_rng);
        cr.pos.y = uniY(m_rng);
        cr.vel   = {0.0f, 0.0f};
    };

    // =========================================================================
    // Main per-carrier loop
    // =========================================================================
    for (auto& c : m_carriers) {

        // ---- Rule 2: depletion-region culling (Painter mode only) ----------
        if (painter) {
            const float u = (c.pos.x - m_topLeft.x) / m_size.x;
            const float v = (c.pos.y - m_topLeft.y) / m_size.y;
            const int i_pre = std::clamp(static_cast<int>(u * W), 0, W - 1);
            const int j_pre = std::clamp(static_cast<int>(v * H), 0, H - 1);
            const double dens_local = c.electron
                ? dd.nDensityAt(i_pre, j_pre, n_i, V_T)
                : dd.pDensityAt(i_pre, j_pre, n_i, V_T);
            if (dens_local < DEPLETION_THRESHOLD) {
                if (teleportToBulk(c)) continue;
                // Painter has no bulk for this species -- fall through to
                // the normal transport; the boundary loop will eventually
                // catch the particle and re-inject it.
            }
        }

        // ---- Cell index after any culling ----------------------------------
        const float u = (c.pos.x - m_topLeft.x) / m_size.x;
        const float v = (c.pos.y - m_topLeft.y) / m_size.y;
        const int i = std::clamp(static_cast<int>(u * W), 0, W - 1);
        const int j = std::clamp(static_cast<int>(v * H), 0, H - 1);

        // ---- Real E-field (V/cm) -------------------------------------------
        const float Ex   = dd.electricFieldX(i, j);
        const float Ey   = dd.electricFieldY(i, j);
        const float Emag = std::sqrt(Ex * Ex + Ey * Ey);

        // ---- Local mu (Matthiessen + Caughey-Thomas; v_sat baked in) -------
        const double N_local = static_cast<double>(dd.donorAt(i, j))
                             + static_cast<double>(dd.acceptorAt(i, j));
        const double mu = c.electron
            ? PhysicsEngine::localMobilityElectron(mat, T, N_local, Emag)
            : PhysicsEngine::localMobilityHole    (mat, T, N_local, Emag);

        // ---- Drift velocity (Sze 1.32): electrons against E, holes with E.
        const float drift_sign = c.electron ? -1.0f : +1.0f;
        const float vx_cms = drift_sign * static_cast<float>(mu) * Ex;
        const float vy_cms = drift_sign * static_cast<float>(mu) * Ey;
        float vx = vx_cms * px_per_cm * tau_vis;
        float vy = vy_cms * px_per_cm * tau_vis;

        // ---- Lorentz / Hall rotation of the drift vector -------------------
        if (magnetic) {
            const float sgn = c.electron ? +1.f : -1.f;
            const float vrx = vx * cs - sgn * vy * sn;
            const float vry = sgn * vx * sn + vy * cs;
            vx = vrx;
            vy = vry;
        }

        // ---- Field-free fallback: small isotropic kick + persistence -------
        if (Emag < E_FIELD_THRESHOLD) {
            vx += tinyKick(m_rng) * dt;
            vy += tinyKick(m_rng) * dt;
            c.vel.x = 0.96f * c.vel.x + 0.04f * vx;
            c.vel.y = 0.96f * c.vel.y + 0.04f * vy;
        } else {
            // High-field zone: drift dominates -- overwrite velocity.
            c.vel.x = vx;
            c.vel.y = vy;
        }

        // ---- Rule 3: damped Einstein diffusion (D = mu V_T) ----------------
        const float D_cm2_per_s = static_cast<float>(mu * V_T);
        const float sigma_cm    = std::sqrt(
            2.0f * std::max(D_cm2_per_s, 0.0f) * std::max(dt_phys, 0.0f));
        const float sigma_px    = sigma_cm * px_per_cm * VISUAL_DIFFUSION_DAMPER;
        const float dx_thermal  = sigma_px * thermal(m_rng);
        const float dy_thermal  = sigma_px * thermal(m_rng);

        // ---- Integrate position --------------------------------------------
        c.pos.x += c.vel.x * dt + dx_thermal;
        c.pos.y += c.vel.y * dt + dy_thermal;

        // ---- Rule 1: absorbing + injecting X-boundaries --------------------
        // Crossing either contact wall re-enters from the OPPOSITE edge,
        // creating the steady-state current illusion: the device is
        // implicitly closed into a circuit. Y-axis still uses a soft clamp
        // because the device geometry is 1D in current flow direction; we
        // don't want carriers wrapping vertically.
        if (c.pos.x < xMin) {
            injectFromOppositeEdge(c, /*exitedRight=*/false);
        } else if (c.pos.x > xMax) {
            injectFromOppositeEdge(c, /*exitedRight=*/true);
        }
        // Soft Y-clamp: stop perpendicular motion at the wall without
        // bouncing (preserves drift direction along x; no jarring rebound).
        if (c.pos.y < yMin) {
            c.pos.y  = yMin;
            c.vel.y  = std::max(0.0f, c.vel.y);
        } else if (c.pos.y > yMax) {
            c.pos.y  = yMax;
            c.vel.y  = std::min(0.0f, c.vel.y);
        }
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
