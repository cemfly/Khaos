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
// Local formatting helpers
// =============================================================================
namespace {

    [[nodiscard]] std::string fix(double v, int prec = 3) {
        std::ostringstream oss;
        oss << std::fixed << std::setprecision(prec) << v;
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


// =============================================================================
// Rendering
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
        "Energy Band Diagram (" + std::string(physics.getMaterial().name) + ")",
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
        axisBox.setFillColor(sf::Color(25, 30, 40));
        axisBox.setOutlineColor(sf::Color(100, 110, 130));
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
