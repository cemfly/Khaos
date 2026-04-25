#include "BandDiagram.hpp"

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

    [[nodiscard]] std::string sci(double v, int prec = 2) {
        std::ostringstream oss;
        oss << std::scientific << std::setprecision(prec) << v;
        return oss.str();
    }

    [[nodiscard]] std::string fix(double v, int prec = 3) {
        std::ostringstream oss;
        oss << std::fixed << std::setprecision(prec) << v;
        return oss.str();
    }

    void drawDashed(sf::RenderTarget& target,
                    float x1, float x2, float y,
                    sf::Color color,
                    float dash = 8.f,
                    float gap  = 6.f)
    {
        sf::VertexArray lines(sf::PrimitiveType::Lines);
        for (float x = x1; x < x2; x += dash + gap) {
            const float xEnd = std::min(x + dash, x2);
            lines.append(sf::Vertex{{x,    y}, color});
            lines.append(sf::Vertex{{xEnd, y}, color});
        }
        target.draw(lines);
    }

    // Upward photon-transition arrow (VB -> CB).
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
// Construction / coordinate mapping
// =============================================================================
BandDiagram::BandDiagram(const sf::Vector2f& topLeft, const sf::Vector2f& size)
    : m_topLeft(topLeft), m_size(size)
{}

float BandDiagram::energyToY(double E) const {
    const float t = static_cast<float>((E - m_Emin) / (m_Emax - m_Emin));
    return m_topLeft.y + m_size.y * (1.f - t);
}


// =============================================================================
// Rendering
// =============================================================================
void BandDiagram::draw(sf::RenderTarget&    target,
                       const PhysicsEngine& physics,
                       const sf::Font&      font) const
{
    // Background
    sf::RectangleShape bg(m_size);
    bg.setPosition(m_topLeft);
    bg.setFillColor(palette::ViewBg);
    bg.setOutlineColor(palette::PanelEdge);
    bg.setOutlineThickness(1.f);
    target.draw(bg);

    const float bandW  = m_size.x * 0.70f;
    const float bandX1 = m_topLeft.x + 20.f;
    const float bandX2 = m_topLeft.x + bandW - 10.f;
    const float fdX1   = m_topLeft.x + bandW + 10.f;
    const float fdX2   = m_topLeft.x + m_size.x - 15.f;

    target.draw(makeText(font, "Energy Band Diagram", 16, palette::TextLight,
                         {m_topLeft.x + 10.f, m_topLeft.y + 6.f}));

    // Pull physics state
    const double     Ev = physics.getValenceBandEdge();
    const double     Ec = physics.getConductionBandEdge();
    const double     Ef = physics.getFermiLevel();
    const double     Ed = physics.getDonorLevel();
    const double     Ea = physics.getAcceptorLevel();
    const DopingType dt = physics.getDopingType();

    // Valence band shading
    {
        sf::RectangleShape r;
        const float yTop = energyToY(Ev);
        const float yBot = energyToY(m_Emin);
        r.setPosition({bandX1, yTop});
        r.setSize({bandX2 - bandX1, yBot - yTop});
        r.setFillColor(sf::Color(40, 100, 160, 120));
        target.draw(r);
    }
    // Conduction band shading
    {
        sf::RectangleShape r;
        const float yTop = energyToY(m_Emax);
        const float yBot = energyToY(Ec);
        r.setPosition({bandX1, yTop});
        r.setSize({bandX2 - bandX1, yBot - yTop});
        r.setFillColor(sf::Color(160, 80, 40, 120));
        target.draw(r);
    }

    // Band edges
    auto drawEdge = [&](double E, sf::Color col,
                        const std::string& label, const std::string& value)
    {
        const float y = energyToY(E);
        sf::VertexArray line(sf::PrimitiveType::Lines);
        line.append(sf::Vertex{{bandX1, y}, col});
        line.append(sf::Vertex{{bandX2, y}, col});
        target.draw(line);
        target.draw(makeText(font, label, 13, col,
                             {bandX1 - 5.f, y - 18.f}));
        target.draw(makeText(font, value, 11, palette::TextDim,
                             {bandX2 - 90.f, y - 16.f}));
    };
    drawEdge(Ec, palette::ConductionBand, "E_c", fix(Ec, 3) + " eV");
    drawEdge(Ev, palette::ValenceBand,    "E_v", fix(Ev, 3) + " eV");

    // Fermi level (dashed)
    {
        const float y = energyToY(Ef);
        drawDashed(target, bandX1, bandX2, y, palette::FermiLine, 9.f, 5.f);
        target.draw(makeText(font, "E_f", 13, palette::FermiLine,
                             {bandX1 - 5.f, y - 18.f}));
        target.draw(makeText(font, fix(Ef, 3) + " eV", 11,
                             sf::Color(220, 220, 160),
                             {bandX2 - 90.f, y - 16.f}));
    }

    // Donor / acceptor level
    if (dt == DopingType::NType) {
        const float y = energyToY(Ed);
        drawDashed(target, bandX1 + 30.f, bandX2 - 30.f, y,
                   palette::DonorLevel, 4.f, 3.f);
        target.draw(makeText(font, "E_d (donor)", 11,
                             sf::Color(120, 240, 160),
                             {bandX1 + 30.f, y - 16.f}));
    }
    else if (dt == DopingType::PType) {
        const float y = energyToY(Ea);
        drawDashed(target, bandX1 + 30.f, bandX2 - 30.f, y,
                   palette::AcceptorLevel, 4.f, 3.f);
        target.draw(makeText(font, "E_a (acceptor)", 11,
                             sf::Color(240, 140, 140),
                             {bandX1 + 30.f, y - 16.f}));
    }

    // ---- Photon absorption arrow (only when hv > E_g) ---------------------
    if (physics.isOpticallyPumped()) {
        const double hv    = physics.getPhotonEnergy();
        const double E_top = std::min<double>(Ev + hv, m_Emax - 0.02);

        const float  x  = bandX1 + 0.22f * (bandX2 - bandX1);
        const float  y0 = energyToY(Ev);
        const float  y1 = energyToY(E_top);
        const sf::Color col = wavelengthColor(physics.getWavelengthNm());

        drawUpwardArrow(target, x, y0, y1, col);
        target.draw(makeText(font, "hv = " + fix(hv, 3) + " eV", 11, col,
                             {x + 10.f, y1 - 4.f}));
    }

    // Energy axis ticks
    {
        sf::VertexArray ticks(sf::PrimitiveType::Lines);
        const sf::Color axisCol(120, 130, 150);
        for (float E = 0.f; E <= m_Emax; E += 0.2f) {
            const float y = energyToY(E);
            ticks.append(sf::Vertex{{bandX1 - 4.f, y}, axisCol});
            ticks.append(sf::Vertex{{bandX1,       y}, axisCol});
        }
        target.draw(ticks);
    }

    // ---- Fermi-Dirac curve panel -----------------------------------------
    target.draw(makeText(font, "f(E)", 13, palette::TextLight,
                         {fdX1, m_topLeft.y + 26.f}));
    {
        sf::RectangleShape axisBox;
        axisBox.setPosition({fdX1, m_topLeft.y + 40.f});
        axisBox.setSize({fdX2 - fdX1, m_size.y - 60.f});
        axisBox.setFillColor(sf::Color(25, 30, 40));
        axisBox.setOutlineColor(sf::Color(100, 110, 130));
        axisBox.setOutlineThickness(1.f);
        target.draw(axisBox);
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
        target.draw(curve);
    }
    {
        const float y = energyToY(Ef);
        const float x = fdX1 + 0.5f * (fdX2 - fdX1);
        sf::CircleShape pt(3.f);
        pt.setOrigin({3.f, 3.f});
        pt.setPosition({x, y});
        pt.setFillColor(palette::FermiLine);
        target.draw(pt);
    }

    // ---- Numeric readouts at bottom of the band panel ---------------------
    const float textX = m_topLeft.x + 10.f;
    float       textY = m_topLeft.y + m_size.y - 110.f;

    auto statusLine = [&](const std::string& s, sf::Color col) {
        target.draw(makeText(font, s, 12, col, {textX, textY}));
        textY += 16.f;
    };

    statusLine("T     = " + fix(physics.getTemperature(), 1) + " K",
               palette::TextDim);
    statusLine("E_g   = " + fix(physics.getBandgap(), 4) + " eV",
               palette::TextDim);
    statusLine("n_i   = " + sci(physics.getIntrinsicCarrier()) + " cm^-3",
               palette::TextDim);
    statusLine("n     = " + sci(physics.getTotalElectronConc()) + " cm^-3",
               palette::Electron);
    statusLine("p     = " + sci(physics.getTotalHoleConc()) + " cm^-3",
               palette::Hole);
    statusLine("sigma = " + sci(physics.getConductivity()) + " S/cm",
               palette::Conductivity);
}
