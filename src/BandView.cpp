// =============================================================================
// BandView.cpp
//
//   Author : dex / cemfly-april2026
//   License: MIT
// =============================================================================

#include "BandView.hpp"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <sstream>
#include <string>

#include "Palette.hpp"

using palette::makeText;
using palette::wavelengthColor;


// =============================================================================
// Local helpers
// =============================================================================
namespace {

[[nodiscard]] std::string fix(double v, int prec = 3) {
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(prec) << v;
    return oss.str();
}

[[nodiscard]] std::string sci(double v, int prec = 2) {
    std::ostringstream oss;
    oss << std::scientific << std::setprecision(prec) << v;
    return oss.str();
}

void drawDashed(sf::RenderTarget& target,
                float x1, float x2, float y,
                sf::Color color,
                float dash = 8.f, float gap = 6.f)
{
    sf::VertexArray lines(sf::PrimitiveType::Lines);
    for (float x = x1; x < x2; x += dash + gap) {
        const float xEnd = std::min(x + dash, x2);
        lines.append(sf::Vertex{{x,    y}, color});
        lines.append(sf::Vertex{{xEnd, y}, color});
    }
    target.draw(lines);
}

void drawUpwardArrow(sf::RenderTarget& target,
                     float x, float yBottom, float yTop,
                     sf::Color color,
                     float headSize = 7.f)
{
    sf::VertexArray shaft(sf::PrimitiveType::Lines);
    shaft.append(sf::Vertex{{x, yBottom}, color});
    shaft.append(sf::Vertex{{x, yTop   }, color});
    target.draw(shaft);

    sf::ConvexShape head;
    head.setPointCount(3);
    head.setPoint(0, {x,             yTop});
    head.setPoint(1, {x - headSize,  yTop + headSize});
    head.setPoint(2, {x + headSize,  yTop + headSize});
    head.setFillColor(color);
    target.draw(head);
}

// Polyline plot helper -- maps a unitless y-buffer into a vertical pixel
// window [yPxTop, yPxBot] using the supplied [yMin, yMax] data range.
// Zero-allocation: takes a pre-allocated SFML VertexArray by reference.
void plotPolyline(sf::RenderTarget&     target,
                  std::span<const float> y_data,
                  float xL, float xR,
                  float yPxTop, float yPxBot,
                  float yMin,   float yMax,
                  sf::Color col, float thickness = 1.0f)
{
    if (y_data.empty() || yMax <= yMin) return;
    const float xspan = xR - xL;
    const float yspan = yPxBot - yPxTop;
    const float invY  = 1.0f / (yMax - yMin);

    sf::VertexArray strip(sf::PrimitiveType::LineStrip);
    const std::size_t N = y_data.size();
    for (std::size_t i = 0; i < N; ++i) {
        const float t  = (N == 1) ? 0.0f
                                  : static_cast<float>(i) /
                                    static_cast<float>(N - 1);
        const float x  = xL + t * xspan;
        const float ny = std::clamp((y_data[i] - yMin) * invY, 0.0f, 1.0f);
        const float y  = yPxBot - ny * yspan;
        strip.append(sf::Vertex{{x, y}, col});
    }
    target.draw(strip);
    // Thick line trick: draw a couple of 1px offset copies. SFML 3 has no
    // native "line thickness" so this is the cheapest way to fatten lines
    // without spawning RectangleShapes per-segment.
    if (thickness > 1.0f) {
        for (int dx = 1; dx < static_cast<int>(thickness); ++dx) {
            sf::VertexArray copy = strip;
            for (std::size_t k = 0; k < copy.getVertexCount(); ++k) {
                copy[k].position.y += static_cast<float>(dx);
            }
            target.draw(copy);
        }
    }
}

} // namespace


// =============================================================================
// Construction
// =============================================================================
BandView::BandView(unsigned int width, unsigned int height) {
    if (!m_rt.resize({width, height})) { /* unrecoverable; check texture size */ }
    m_rt.setSmooth(true);
}

float BandView::energyToY(double E) const noexcept {
    const float t = static_cast<float>((E - m_Emin) / (m_Emax - m_Emin));
    return static_cast<float>(m_rt.getSize().y) * (1.f - t);
}

float BandView::xCellToPx(int i, int W_cells,
                          float xL, float xR) const noexcept
{
    if (W_cells <= 1) return xL;
    const float t = static_cast<float>(i) / static_cast<float>(W_cells - 1);
    return xL + t * (xR - xL);
}

void BandView::ensureBuffersFor(int W_cells) {
    if (W_cells <= 0) return;
    const auto need = static_cast<std::size_t>(W_cells);
    if (m_psi_cut.size()    < need) m_psi_cut   .resize(need);
    if (m_Ec_cut.size()     < need) m_Ec_cut    .resize(need);
    if (m_Ev_cut.size()     < need) m_Ev_cut    .resize(need);
    if (m_Efield_cut.size() < need) m_Efield_cut.resize(need);
    if (m_doping_cut.size() < need) m_doping_cut.resize(need);
    if (m_phiN_cut.size()   < need) m_phiN_cut  .resize(need);
    if (m_phiP_cut.size()   < need) m_phiP_cut  .resize(need);
}


// =============================================================================
// Flat rendering (legacy)
// =============================================================================
void BandView::render(const PhysicsEngine& physics, const sf::Font& font) {
    m_rt.clear(palette::ViewBg);

    const float W = static_cast<float>(m_rt.getSize().x);
    const float H = static_cast<float>(m_rt.getSize().y);

    const float bandW  = W * 0.70f;
    const float bandX1 = 20.f;
    const float bandX2 = bandW - 10.f;
    const float fdX1   = bandW + 10.f;
    const float fdX2   = W - 15.f;

    m_rt.draw(makeText(font,
        "Energy Band Diagram (" + std::string(physics.getMaterial().name)
        + ") -- flat",
        16, palette::TextLight, {10.f, 6.f}));

    const double     Ev = physics.getValenceBandEdge();
    const double     Ec = physics.getConductionBandEdge();
    const double     Ef = physics.getFermiLevel();
    const double     Ed = physics.getDonorLevel();
    const double     Ea = physics.getAcceptorLevel();
    const DopingType dt = physics.getDopingType();

    // ---- Band shading -----------------------------------------------------
    {
        sf::RectangleShape r;
        const float yTop = energyToY(Ev);
        const float yBot = energyToY(m_Emin);
        r.setPosition({bandX1, yTop});
        r.setSize({bandX2 - bandX1, yBot - yTop});
        r.setFillColor(sf::Color(40, 100, 160, 120));
        m_rt.draw(r);
    }
    {
        sf::RectangleShape r;
        const float yTop = energyToY(m_Emax);
        const float yBot = energyToY(Ec);
        r.setPosition({bandX1, yTop});
        r.setSize({bandX2 - bandX1, yBot - yTop});
        r.setFillColor(sf::Color(160, 80, 40, 120));
        m_rt.draw(r);
    }

    // ---- Band edges -------------------------------------------------------
    auto drawEdge = [&](double E, sf::Color col,
                        const std::string& label, const std::string& value)
    {
        const float y = energyToY(E);
        sf::VertexArray line(sf::PrimitiveType::Lines);
        line.append(sf::Vertex{{bandX1, y}, col});
        line.append(sf::Vertex{{bandX2, y}, col});
        m_rt.draw(line);
        m_rt.draw(makeText(font, label, 13, col,
                           {bandX1 - 5.f, y - 18.f}));
        m_rt.draw(makeText(font, value, 11, palette::TextDim,
                           {bandX2 - 90.f, y - 16.f}));
    };
    drawEdge(Ec, palette::ConductionBand, "E_c", fix(Ec, 3) + " eV");
    drawEdge(Ev, palette::ValenceBand,    "E_v", fix(Ev, 3) + " eV");

    // ---- Fermi level (dashed) ---------------------------------------------
    {
        const float y = energyToY(Ef);
        drawDashed(m_rt, bandX1, bandX2, y, palette::FermiLine, 9.f, 5.f);
        m_rt.draw(makeText(font, "E_f", 13, palette::FermiLine,
                           {bandX1 - 5.f, y - 18.f}));
        m_rt.draw(makeText(font, fix(Ef, 3) + " eV", 11,
                           sf::Color(220, 220, 160),
                           {bandX2 - 90.f, y - 16.f}));
    }

    // ---- Donor / acceptor level -------------------------------------------
    if (dt == DopingType::NType) {
        const float y = energyToY(Ed);
        drawDashed(m_rt, bandX1 + 30.f, bandX2 - 30.f, y,
                   palette::DonorLevel, 4.f, 3.f);
        m_rt.draw(makeText(font,
            "E_d (" + std::string(physics.getMaterial().donorSpecies) + ")",
            11, sf::Color(120, 240, 160),
            {bandX1 + 30.f, y - 16.f}));
    } else if (dt == DopingType::PType) {
        const float y = energyToY(Ea);
        drawDashed(m_rt, bandX1 + 30.f, bandX2 - 30.f, y,
                   palette::AcceptorLevel, 4.f, 3.f);
        m_rt.draw(makeText(font,
            "E_a (" + std::string(physics.getMaterial().acceptorSpecies) + ")",
            11, sf::Color(240, 140, 140),
            {bandX1 + 30.f, y - 16.f}));
    }

    // ---- Photon arrow (when hv > E_g) -------------------------------------
    if (physics.isOpticallyPumped()) {
        const double hv    = physics.getPhotonEnergy();
        const double E_top = std::min<double>(Ev + hv, m_Emax - 0.02);

        const float  x  = bandX1 + 0.22f * (bandX2 - bandX1);
        const float  y0 = energyToY(Ev);
        const float  y1 = energyToY(E_top);
        const sf::Color col = wavelengthColor(physics.getWavelengthNm());

        drawUpwardArrow(m_rt, x, y0, y1, col);
        m_rt.draw(makeText(font, "hv = " + fix(hv, 3) + " eV", 11, col,
                           {x + 10.f, y1 - 4.f}));
    }

    // ---- Energy axis ticks ------------------------------------------------
    {
        sf::VertexArray ticks(sf::PrimitiveType::Lines);
        const sf::Color axisCol(120, 130, 150);
        for (float E = 0.f; E <= m_Emax; E += 0.2f) {
            const float y = energyToY(E);
            ticks.append(sf::Vertex{{bandX1 - 4.f, y}, axisCol});
            ticks.append(sf::Vertex{{bandX1,       y}, axisCol});
        }
        m_rt.draw(ticks);
    }

    // ---- Fermi-Dirac panel ------------------------------------------------
    m_rt.draw(makeText(font, "f(E)", 13, palette::TextLight,
                       {fdX1, 26.f}));
    {
        sf::RectangleShape axisBox;
        axisBox.setPosition({fdX1, 40.f});
        axisBox.setSize({fdX2 - fdX1, H - 60.f});
        axisBox.setFillColor(sf::Color(255, 255, 255));
        axisBox.setOutlineColor(palette::PanelEdge);
        axisBox.setOutlineThickness(1.f);
        m_rt.draw(axisBox);
    }
    {
        sf::VertexArray curve(sf::PrimitiveType::LineStrip);
        const int samples = 120;
        for (int i = 0; i <= samples; ++i) {
            const double E = m_Emin + (m_Emax - m_Emin) * i / samples;
            const double f = physics.fermiDirac(E);
            const float  y = energyToY(E);
            const float  x = fdX1 + static_cast<float>(f) * (fdX2 - fdX1);
            curve.append(sf::Vertex{{x, y}, palette::FermiLine});
        }
        m_rt.draw(curve);
    }
    {
        const float y = energyToY(Ef);
        const float x = fdX1 + 0.5f * (fdX2 - fdX1);
        sf::CircleShape pt(3.f);
        pt.setOrigin({3.f, 3.f});
        pt.setPosition({x, y});
        pt.setFillColor(palette::FermiLine);
        m_rt.draw(pt);
    }

    m_rt.display();
}


// =============================================================================
// Spatial rendering [Phase 2]
// -----------------------------------------------------------------------------
// Pipeline:
//   1. Pull Ec(x), Ev(x), |E|(x), Nd-Na(x) along the chosen cut row from
//      the drift-diffusion grid (zero-allocation: pre-sized internal
//      buffers).
//   2. Auto-fit the energy axis to (min Ev - 0.2, max Ec + 0.2) so the
//      bent bands stay on screen even at large built-in potentials.
//   3. Paint the doping background tint, draw Ec/Ev/Ef polylines, then
//      overlay the |E|(x) curve in a smaller secondary panel below the
//      band canvas.
// =============================================================================
void BandView::renderSpatial(const PhysicsEngine&  physics,
                             const DriftDiffusion& dd,
                             int                   j_cut,
                             const sf::Font&       font)
{
    m_rt.clear(palette::ViewBg);

    const float W = static_cast<float>(m_rt.getSize().x);
    const float H = static_cast<float>(m_rt.getSize().y);

    // Layout: top 70% = bands, bottom 30% = E-field overlay.
    const float xL = 50.f;
    const float xR = W - 25.f;
    const float yBandsTop = 40.f;
    const float yBandsBot = H * 0.66f;
    const float yEfTop    = H * 0.72f;
    const float yEfBot    = H - 30.f;

    const int W_cells = dd.width();
    if (W_cells <= 0) {
        m_rt.draw(makeText(font, "(no grid)", 16, palette::TextLight,
                           {10.f, 10.f}));
        m_rt.display();
        return;
    }
    ensureBuffersFor(W_cells);

    j_cut = std::clamp(j_cut, 0, dd.height() - 1);

    // Reference Ec0/Ev0 from the bulk physics (the painter Painter mode
    // does not have a global Ef; we still draw the engine's Fermi level
    // as a constant reference line for orientation).
    const double Ev0 = physics.getValenceBandEdge();      // 0
    const double Ec0 = physics.getConductionBandEdge();   // = Eg
    const double Ef  = physics.getFermiLevel();

    // Pull cuts -- zero-allocation thanks to ensureBuffersFor() above.
    dd.sampleHorizontalCut(j_cut,
        std::span<float>(m_psi_cut.data(), W_cells),
        std::span<float>(m_Ec_cut .data(), W_cells),
        std::span<float>(m_Ev_cut .data(), W_cells),
        Ec0, Ev0);
    dd.sampleEFieldCut(j_cut,
        std::span<float>(m_Efield_cut.data(), W_cells));
    dd.sampleDopingCut(j_cut,
        std::span<float>(m_doping_cut.data(), W_cells));
    dd.sampleQuasiFermiCut(j_cut,
        std::span<float>(m_phiN_cut.data(), W_cells),
        std::span<float>(m_phiP_cut.data(), W_cells));

    // ---- Auto-fit the energy axis to the bent profile -------------------
    float Emin =  static_cast<float>(Ev0 - 0.20);
    float Emax =  static_cast<float>(Ec0 + 0.20);
    for (int i = 0; i < W_cells; ++i) {
        Emin = std::min(Emin, m_Ev_cut[i] - 0.05f);
        Emax = std::max(Emax, m_Ec_cut[i] + 0.05f);
    }
    m_Emin = Emin;
    m_Emax = Emax;

    // Map to band-region pixels (override the default texture-height map).
    auto eToY = [&](float E) {
        const float t = (E - Emin) / (Emax - Emin);
        return yBandsBot - t * (yBandsBot - yBandsTop);
    };

    // ---- Title + cut info -----------------------------------------------
    m_rt.draw(makeText(font,
        "Spatial Band Diagram (" + std::string(physics.getMaterial().name)
        + ")  --  cut row j = " + std::to_string(j_cut),
        16, palette::TextLight, {10.f, 6.f}));

    m_rt.draw(makeText(font,
        "x: 0 -> " + fix(static_cast<double>(W_cells)
                         * dd.cellPitchCm() * 1.0e4, 1) + " um",
        11, palette::TextDim, {10.f, 24.f}));

    // ---- Doping background tint -----------------------------------------
    // Blue stripe under N-doped cells, orange under P-doped cells. Helps
    // the viewer locate the metallurgical junction.
    {
        for (int i = 0; i < W_cells; ++i) {
            const float net = m_doping_cut[i];
            if (net == 0.0f) continue;
            const float lg = std::log10(std::abs(net) + 1.0f) / 20.0f;
            const auto a = static_cast<std::uint8_t>(
                std::clamp(lg * 80.0f, 0.0f, 80.0f));
            sf::RectangleShape r;
            const float x0 = xCellToPx(i,     W_cells, xL, xR);
            const float x1 = xCellToPx(i + 1, W_cells, xL, xR);
            r.setPosition({x0, yBandsTop});
            r.setSize({std::max(1.0f, x1 - x0), yBandsBot - yBandsTop});
            r.setFillColor(net > 0 ? sf::Color(60, 140, 255, a)
                                   : sf::Color(255, 140, 60, a));
            m_rt.draw(r);
        }
    }

    // ---- Energy axis grid -----------------------------------------------
    {
        const sf::Color axisCol(120, 130, 150, 140);
        sf::VertexArray ticks(sf::PrimitiveType::Lines);
        for (float E = std::ceil(Emin * 5.0f) / 5.0f; E <= Emax; E += 0.2f) {
            const float y = eToY(E);
            ticks.append(sf::Vertex{{xL - 6.f, y}, axisCol});
            ticks.append(sf::Vertex{{xR,        y}, axisCol});
            m_rt.draw(makeText(font, fix(E, 1), 9, palette::TextDim,
                               {xL - 38.f, y - 8.f}));
        }
        m_rt.draw(ticks);
    }

    // ---- Quasi-Fermi levels [Phase 3] -----------------------------------
    //
    //   E_fn(x) = E_i_ref - phi_n(x)
    //   E_fp(x) = E_i_ref - phi_p(x)
    //
    // In painter mode there is no globally-defined Fermi level (doping is
    // spatial), so we use the intrinsic-level reference E_i_ref = Eg/2.
    // In equilibrium phi_n = phi_p = 0 -> both lines collapse to E_i_ref
    // and the diagram reproduces the classical equilibrium picture. Under
    // bias the lines split inside the depletion region by ~ q V_a (Sze
    // Fig. 2.10 / Pierret Fig. 5.10).
    const float E_i_ref = static_cast<float>(0.5 * (Ec0 + Ev0));
    const bool  has_bias = std::abs(dd.appliedBias()) > 1.0e-4f;

    if (has_bias) {
        // E_fn (blue dashed) -- electron quasi-Fermi level
        sf::Color  ef_n_col(110, 175, 255);
        {
            sf::VertexArray strip(sf::PrimitiveType::LineStrip);
            for (int i = 0; i < W_cells; ++i) {
                const float Efn = E_i_ref - m_phiN_cut[i];
                strip.append(sf::Vertex{
                    {xCellToPx(i, W_cells, xL, xR), eToY(Efn)}, ef_n_col});
            }
            // Dashed feel: draw every other segment by varying alpha; the
            // simplest way is to plot the polyline directly and let the
            // colour distinguish it from the solid Ec/Ev curves.
            m_rt.draw(strip);
            sf::VertexArray strip2 = strip;
            for (std::size_t k = 0; k < strip2.getVertexCount(); ++k)
                strip2[k].position.y += 1.0f;
            m_rt.draw(strip2);

            const float y_first = eToY(E_i_ref - m_phiN_cut.front());
            m_rt.draw(makeText(font, "E_fn", 11, ef_n_col,
                               {xL - 38.f, y_first - 16.f}));
        }
        // E_fp (orange dashed) -- hole quasi-Fermi level
        sf::Color ef_p_col(255, 165, 90);
        {
            sf::VertexArray strip(sf::PrimitiveType::LineStrip);
            for (int i = 0; i < W_cells; ++i) {
                const float Efp = E_i_ref - m_phiP_cut[i];
                strip.append(sf::Vertex{
                    {xCellToPx(i, W_cells, xL, xR), eToY(Efp)}, ef_p_col});
            }
            m_rt.draw(strip);
            sf::VertexArray strip2 = strip;
            for (std::size_t k = 0; k < strip2.getVertexCount(); ++k)
                strip2[k].position.y += 1.0f;
            m_rt.draw(strip2);

            const float y_first = eToY(E_i_ref - m_phiP_cut.front());
            m_rt.draw(makeText(font, "E_fp", 11, ef_p_col,
                               {xL - 38.f, y_first + 4.f}));
        }
        // Bias readout in the corner.
        m_rt.draw(makeText(font,
            "V_a = " + fix(static_cast<double>(dd.appliedBias()), 3) + " V",
            12, sf::Color(255, 230, 130),
            {xR - 100.f, yBandsTop + 20.f}));

        // Highlight the split E_fn - E_fp at the painted junction (sign
        // change in psi). Draw a vertical bracket showing the split.
        for (int i = 1; i < W_cells; ++i) {
            const float a = m_psi_cut[i];
            const float b = m_psi_cut[i - 1];
            if ((a > 0) != (b > 0)) {
                const float x  = xCellToPx(i, W_cells, xL, xR);
                const float y1 = eToY(E_i_ref - m_phiN_cut[i]);
                const float y2 = eToY(E_i_ref - m_phiP_cut[i]);
                sf::VertexArray bracket(sf::PrimitiveType::Lines);
                bracket.append(sf::Vertex{{x, y1},
                    sf::Color(255, 240, 120, 200)});
                bracket.append(sf::Vertex{{x, y2},
                    sf::Color(255, 240, 120, 200)});
                m_rt.draw(bracket);
                const double split_eV = static_cast<double>(
                    m_phiP_cut[i] - m_phiN_cut[i]);
                m_rt.draw(makeText(font,
                    "  qV = " + fix(split_eV, 3) + " eV", 10,
                    sf::Color(255, 240, 120),
                    {x + 4.f, std::min(y1, y2) - 6.f}));
                break;   // one annotation is enough
            }
        }
    } else {
        // Equilibrium fallback: the engine's global Ef as a single line.
        const float y = eToY(static_cast<float>(Ef));
        drawDashed(m_rt, xL, xR, y, palette::FermiLine, 9.f, 5.f);
        m_rt.draw(makeText(font, "E_f", 13, palette::FermiLine,
                           {xR + 4.f, y - 8.f}));
    }

    // ---- Bent E_c / E_v polylines ---------------------------------------
    {
        sf::VertexArray ec(sf::PrimitiveType::LineStrip);
        sf::VertexArray ev(sf::PrimitiveType::LineStrip);
        for (int i = 0; i < W_cells; ++i) {
            const float x = xCellToPx(i, W_cells, xL, xR);
            ec.append(sf::Vertex{{x, eToY(m_Ec_cut[i])},
                                 palette::ConductionBand});
            ev.append(sf::Vertex{{x, eToY(m_Ev_cut[i])},
                                 palette::ValenceBand});
        }
        // Light "shadow" copies for thicker line feel.
        m_rt.draw(ec);
        m_rt.draw(ev);
        for (int dy = 1; dy <= 1; ++dy) {
            sf::VertexArray ec2 = ec, ev2 = ev;
            for (std::size_t k = 0; k < ec2.getVertexCount(); ++k)
                ec2[k].position.y += static_cast<float>(dy);
            for (std::size_t k = 0; k < ev2.getVertexCount(); ++k)
                ev2[k].position.y += static_cast<float>(dy);
            m_rt.draw(ec2);
            m_rt.draw(ev2);
        }
    }

    // ---- Edge labels ----------------------------------------------------
    {
        const float y_ec_left = eToY(m_Ec_cut.front());
        const float y_ev_left = eToY(m_Ev_cut.front());
        m_rt.draw(makeText(font, "E_c", 13, palette::ConductionBand,
                           {xL - 28.f, y_ec_left - 8.f}));
        m_rt.draw(makeText(font, "E_v", 13, palette::ValenceBand,
                           {xL - 28.f, y_ev_left - 8.f}));

        const float y_ec_right = eToY(m_Ec_cut.back());
        const float y_ev_right = eToY(m_Ev_cut.back());
        const double V_bi = static_cast<double>(
            m_Ec_cut.front() - m_Ec_cut.back());
        m_rt.draw(makeText(font,
            "V_bi ~ " + fix(V_bi, 3) + " V", 11,
            sf::Color(220, 230, 200),
            {xR - 110.f, y_ec_right - 18.f}));
        (void)y_ev_right;
    }

    // ---- E-field secondary panel ----------------------------------------
    {
        sf::RectangleShape box;
        box.setPosition({xL, yEfTop});
        box.setSize({xR - xL, yEfBot - yEfTop});
        box.setFillColor(sf::Color(15, 18, 24, 200));
        box.setOutlineColor(palette::PanelEdge);
        box.setOutlineThickness(1.f);
        m_rt.draw(box);

        m_rt.draw(makeText(font, "|E|(x)  [V/cm]", 11, palette::TextLight,
                           {xL + 6.f, yEfTop + 4.f}));

        // Scale: peak of cut, padded.
        float Emax_field = 1.0f;
        for (int i = 0; i < W_cells; ++i)
            Emax_field = std::max(Emax_field, m_Efield_cut[i]);
        const float Epeak_padded = Emax_field * 1.1f;

        plotPolyline(m_rt,
            std::span<const float>(m_Efield_cut.data(), W_cells),
            xL + 4.f, xR - 4.f,
            yEfTop + 18.f, yEfBot - 6.f,
            0.0f, Epeak_padded,
            sf::Color(255, 200, 80), 2.0f);

        m_rt.draw(makeText(font, "peak: " + sci(Emax_field, 2)
                          + " V/cm", 10, palette::TextDim,
                          {xR - 130.f, yEfTop + 4.f}));
    }

    m_rt.display();
}
