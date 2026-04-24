#include "BandDiagram.hpp"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <sstream>
#include <string>


// =============================================================================
// Local helpers
// =============================================================================
namespace {

    std::string formatSci(double v, int precision = 2) {
        std::ostringstream oss;
        oss << std::scientific << std::setprecision(precision) << v;
        return oss.str();
    }

    std::string formatFixed(double v, int precision = 3) {
        std::ostringstream oss;
        oss << std::fixed << std::setprecision(precision) << v;
        return oss.str();
    }

    sf::Text makeLabel(const std::string& str,
                       const sf::Font&    font,
                       unsigned int       size,
                       sf::Color          color,
                       const sf::Vector2f& pos)
    {
        sf::Text t;
        t.setFont(font);
        t.setString(str);
        t.setCharacterSize(size);
        t.setFillColor(color);
        t.setPosition(pos);
        return t;
    }

    void drawDashed(sf::RenderTarget& target,
                    float x1, float x2, float y,
                    sf::Color color,
                    float dash = 8.f,
                    float gap  = 6.f)
    {
        sf::VertexArray lines(sf::Lines);
        for (float x = x1; x < x2; x += dash + gap) {
            const float xEnd = std::min(x + dash, x2);
            lines.append(sf::Vertex({x,    y}, color));
            lines.append(sf::Vertex({xEnd, y}, color));
        }
        target.draw(lines);
    }

    // Draws a simple upward arrow (photon transition  VB -> CB).
    void drawUpwardArrow(sf::RenderTarget& target,
                         float x, float yBottom, float yTop,
                         sf::Color color,
                         float headSize = 7.f)
    {
        sf::VertexArray shaft(sf::Lines);
        shaft.append(sf::Vertex({x, yBottom}, color));
        shaft.append(sf::Vertex({x, yTop   }, color));
        target.draw(shaft);

        sf::ConvexShape head;
        head.setPointCount(3);
        head.setPoint(0, {x,             yTop});
        head.setPoint(1, {x - headSize,  yTop + headSize});
        head.setPoint(2, {x + headSize,  yTop + headSize});
        head.setFillColor(color);
        target.draw(head);
    }

    // Same wavelength -> RGB mapping used in the UI panel; kept local to
    // avoid an extra header dependency.
    sf::Color wavelengthToColor(double lambda_nm) {
        if (lambda_nm < 380.0 || lambda_nm > 780.0)
            return sf::Color(180, 180, 200);
        if (lambda_nm < 440.0) return sf::Color(130, 100, 255);
        if (lambda_nm < 490.0) return sf::Color( 90, 180, 255);
        if (lambda_nm < 510.0) return sf::Color( 80, 255, 220);
        if (lambda_nm < 580.0) return sf::Color(120, 255, 120);
        if (lambda_nm < 645.0) return sf::Color(255, 220,  80);
        return                        sf::Color(255, 110,  90);
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
void BandDiagram::draw(sf::RenderTarget&   target,
                       const PhysicsEngine& physics,
                       const sf::Font&      font) const
{
    // Background
    sf::RectangleShape bg(m_size);
    bg.setPosition(m_topLeft);
    bg.setFillColor(sf::Color(18, 22, 30));
    bg.setOutlineColor(sf::Color(80, 90, 110));
    bg.setOutlineThickness(1.f);
    target.draw(bg);

    const float bandW  = m_size.x * 0.70f;
    const float bandX1 = m_topLeft.x + 20.f;
    const float bandX2 = m_topLeft.x + bandW - 10.f;
    const float fdX1   = m_topLeft.x + bandW + 10.f;
    const float fdX2   = m_topLeft.x + m_size.x - 15.f;

    target.draw(makeLabel("Energy Band Diagram", font, 16,
                          sf::Color(220, 220, 230),
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
        r.setPosition(bandX1, yTop);
        r.setSize({bandX2 - bandX1, yBot - yTop});
        r.setFillColor(sf::Color(40, 100, 160, 120));
        target.draw(r);
    }
    // Conduction band shading
    {
        sf::RectangleShape r;
        const float yTop = energyToY(m_Emax);
        const float yBot = energyToY(Ec);
        r.setPosition(bandX1, yTop);
        r.setSize({bandX2 - bandX1, yBot - yTop});
        r.setFillColor(sf::Color(160, 80, 40, 120));
        target.draw(r);
    }

    // Band edges
    auto drawEdge = [&](double E, sf::Color col,
                        const std::string& label, const std::string& value)
    {
        const float y = energyToY(E);
        sf::VertexArray line(sf::Lines);
        line.append(sf::Vertex({bandX1, y}, col));
        line.append(sf::Vertex({bandX2, y}, col));
        target.draw(line);
        target.draw(makeLabel(label, font, 13, col,
                              {bandX1 - 5.f, y - 18.f}));
        target.draw(makeLabel(value, font, 11, sf::Color(200, 200, 210),
                              {bandX2 - 90.f, y - 16.f}));
    };
    drawEdge(Ec, sf::Color(230, 140, 100),
             "E_c", formatFixed(Ec, 3) + " eV");
    drawEdge(Ev, sf::Color(120, 180, 230),
             "E_v", formatFixed(Ev, 3) + " eV");

    // Fermi level (dashed)
    {
        const float y = energyToY(Ef);
        drawDashed(target, bandX1, bandX2, y,
                   sf::Color(255, 240, 90), 9.f, 5.f);
        target.draw(makeLabel("E_f", font, 13,
                              sf::Color(255, 240, 90),
                              {bandX1 - 5.f, y - 18.f}));
        target.draw(makeLabel(formatFixed(Ef, 3) + " eV", font, 11,
                              sf::Color(220, 220, 160),
                              {bandX2 - 90.f, y - 16.f}));
    }

    // Donor / acceptor level
    if (dt == DopingType::NType) {
        const float y = energyToY(Ed);
        drawDashed(target, bandX1 + 30.f, bandX2 - 30.f, y,
                   sf::Color(80, 230, 140), 4.f, 3.f);
        target.draw(makeLabel("E_d (donor)", font, 11,
                              sf::Color(120, 240, 160),
                              {bandX1 + 30.f, y - 16.f}));
    }
    else if (dt == DopingType::PType) {
        const float y = energyToY(Ea);
        drawDashed(target, bandX1 + 30.f, bandX2 - 30.f, y,
                   sf::Color(240, 120, 120), 4.f, 3.f);
        target.draw(makeLabel("E_a (acceptor)", font, 11,
                              sf::Color(240, 140, 140),
                              {bandX1 + 30.f, y - 16.f}));
    }

    // ---- Photon absorption arrow (only when hv > E_g) ---------------------
    //
    // Visualises the interband transition  E_v -> E_c  triggered by an
    // incident photon of energy hv. The arrow starts at the valence band
    // edge and climbs by hv (clamped to the plotting range). Its colour
    // matches the selected wavelength.
    //
    // Reference: Kittel Ch. 15, Kasap Ch. 5.14.
    if (physics.isOpticallyPumped()) {
        const double hv = physics.getPhotonEnergy();
        const double E_top = std::min<double>(Ev + hv, m_Emax - 0.02);

        const float  x  = bandX1 + 0.22f * (bandX2 - bandX1);
        const float  y0 = energyToY(Ev);
        const float  y1 = energyToY(E_top);
        const sf::Color col = wavelengthToColor(physics.getWavelengthNm());

        drawUpwardArrow(target, x, y0, y1, col);
        target.draw(makeLabel("hv = " + formatFixed(hv, 3) + " eV",
                              font, 11, col,
                              {x + 10.f, y1 - 4.f}));
    }

    // Energy axis ticks
    {
        sf::VertexArray ticks(sf::Lines);
        const sf::Color axisCol(120, 130, 150);
        for (float E = 0.f; E <= m_Emax; E += 0.2f) {
            const float y = energyToY(E);
            ticks.append(sf::Vertex({bandX1 - 4.f, y}, axisCol));
            ticks.append(sf::Vertex({bandX1,       y}, axisCol));
        }
        target.draw(ticks);
    }

    // ---- Fermi-Dirac curve panel -----------------------------------------
    target.draw(makeLabel("f(E)", font, 13, sf::Color(220, 220, 230),
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
        sf::VertexArray curve(sf::LineStrip);
        const int samples = 120;
        for (int i = 0; i <= samples; ++i) {
            const double E = m_Emin + (m_Emax - m_Emin) * i / samples;
            const double f = physics.fermiDirac(E);
            const float  y = energyToY(E);
            const float  x = fdX1 + static_cast<float>(f) * (fdX2 - fdX1);
            curve.append(sf::Vertex({x, y}, sf::Color(255, 240, 90)));
        }
        target.draw(curve);
    }
    {
        const float y = energyToY(Ef);
        const float x = fdX1 + 0.5f * (fdX2 - fdX1);
        sf::CircleShape pt(3.f);
        pt.setOrigin(3.f, 3.f);
        pt.setPosition(x, y);
        pt.setFillColor(sf::Color(255, 240, 90));
        target.draw(pt);
    }

    // ---- Numeric readouts -------------------------------------------------
    const float textX = m_topLeft.x + 10.f;
    float       textY = m_topLeft.y + m_size.y - 110.f;
    const auto  col   = sf::Color(210, 215, 225);

    target.draw(makeLabel("T     = " + formatFixed(physics.getTemperature(), 1) + " K",
                          font, 12, col, {textX, textY})); textY += 16.f;
    target.draw(makeLabel("E_g   = " + formatFixed(physics.getBandgap(), 4) + " eV",
                          font, 12, col, {textX, textY})); textY += 16.f;
    target.draw(makeLabel("n_i   = " + formatSci(physics.getIntrinsicCarrier()) + " cm^-3",
                          font, 12, col, {textX, textY})); textY += 16.f;
    target.draw(makeLabel("n     = " + formatSci(physics.getTotalElectronConc()) + " cm^-3",
                          font, 12, sf::Color(255, 230, 80), {textX, textY})); textY += 16.f;
    target.draw(makeLabel("p     = " + formatSci(physics.getTotalHoleConc()) + " cm^-3",
                          font, 12, sf::Color(255, 100, 200), {textX, textY})); textY += 16.f;
    target.draw(makeLabel("sigma = " + formatSci(physics.getConductivity()) + " S/cm",
                          font, 12, sf::Color(140, 255, 180), {textX, textY}));
}
