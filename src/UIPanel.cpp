#include "UIPanel.hpp"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <sstream>


// =============================================================================
// Local formatting helpers
// =============================================================================
namespace {

    std::string formatSci(double v, int prec = 2) {
        std::ostringstream oss;
        oss << std::scientific << std::setprecision(prec) << v;
        return oss.str();
    }

    std::string formatFixed(double v, int prec = 2) {
        std::ostringstream oss;
        oss << std::fixed << std::setprecision(prec) << v;
        return oss.str();
    }

    sf::Text makeText(const std::string& s, const sf::Font& f,
                      unsigned size, sf::Color col,
                      const sf::Vector2f& pos)
    {
        sf::Text t;
        t.setFont(f);
        t.setString(s);
        t.setCharacterSize(size);
        t.setFillColor(col);
        t.setPosition(pos);
        return t;
    }

    // Wavelength -> RGB colour (linear pedagogical mapping across the visible
    // spectrum; UV and IR fall back to grey).
    sf::Color wavelengthToColor(double lambda_nm) {
        if (lambda_nm < 380.0 || lambda_nm > 780.0)
            return sf::Color(160, 160, 180);
        if (lambda_nm < 440.0) return sf::Color(130, 100, 255);  // violet
        if (lambda_nm < 490.0) return sf::Color( 90, 180, 255);  // blue
        if (lambda_nm < 510.0) return sf::Color( 80, 255, 220);  // cyan
        if (lambda_nm < 580.0) return sf::Color(120, 255, 120);  // green
        if (lambda_nm < 645.0) return sf::Color(255, 220,  80);  // yellow
        return                        sf::Color(255, 110,  90);  // red
    }

} // namespace


// =============================================================================
// Slider
// =============================================================================
Slider::Slider(const sf::Vector2f& pos,
               float length,
               float minVal,
               float maxVal,
               float initVal,
               bool  logarithmic)
    : m_pos(pos), m_length(length)
    , m_min(minVal), m_max(maxVal)
    , m_value(std::clamp(initVal, minVal, maxVal))
    , m_log(logarithmic)
{}

sf::FloatRect Slider::bounds() const {
    return { m_pos.x - 6.f, m_pos.y - 12.f, m_length + 12.f, 24.f };
}

float Slider::normalized() const {
    if (m_log) {
        const float lmin = std::log10(m_min);
        const float lmax = std::log10(m_max);
        return (std::log10(m_value) - lmin) / (lmax - lmin);
    }
    return (m_value - m_min) / (m_max - m_min);
}

void Slider::setFromNormalized(float t) {
    t = std::clamp(t, 0.f, 1.f);
    float v;
    if (m_log) {
        const float lmin = std::log10(m_min);
        const float lmax = std::log10(m_max);
        v = std::pow(10.f, lmin + t * (lmax - lmin));
    } else {
        v = m_min + t * (m_max - m_min);
    }
    if (std::abs(v - m_value) > 1e-9f) {
        m_value   = v;
        m_changed = true;
    }
}

void Slider::handleEvent(const sf::Event& ev, const sf::Vector2f& mp) {
    if (ev.type == sf::Event::MouseButtonPressed &&
        ev.mouseButton.button == sf::Mouse::Left)
    {
        if (bounds().contains(mp)) {
            m_dragging = true;
            setFromNormalized((mp.x - m_pos.x) / m_length);
        }
    }
    else if (ev.type == sf::Event::MouseButtonReleased &&
             ev.mouseButton.button == sf::Mouse::Left)
    {
        m_dragging = false;
    }
    else if (ev.type == sf::Event::MouseMoved && m_dragging) {
        setFromNormalized((mp.x - m_pos.x) / m_length);
    }
}

void Slider::draw(sf::RenderTarget&  target,
                  const std::string& label,
                  const std::string& valueStr,
                  const sf::Font&    font) const
{
    sf::RectangleShape track({m_length, 4.f});
    track.setPosition(m_pos.x, m_pos.y - 2.f);
    track.setFillColor(sf::Color(90, 100, 120));
    target.draw(track);

    const float handleX = m_pos.x + normalized() * m_length;
    sf::CircleShape handle(8.f);
    handle.setOrigin(8.f, 8.f);
    handle.setPosition(handleX, m_pos.y);
    handle.setFillColor(sf::Color(255, 240, 90));
    handle.setOutlineColor(sf::Color(255, 255, 200));
    handle.setOutlineThickness(1.5f);
    target.draw(handle);

    target.draw(makeText(label, font, 13, sf::Color(210, 215, 225),
                         {m_pos.x, m_pos.y - 28.f}));
    target.draw(makeText(valueStr, font, 13, sf::Color(255, 240, 150),
                         {m_pos.x + m_length + 14.f, m_pos.y - 10.f}));
}


// =============================================================================
// Button
// =============================================================================
Button::Button(const sf::Vector2f& pos,
               const sf::Vector2f& size,
               std::string         label)
    : m_pos(pos), m_size(size), m_label(std::move(label))
{}

void Button::handleEvent(const sf::Event& ev, const sf::Vector2f& mp) {
    if (ev.type == sf::Event::MouseButtonPressed &&
        ev.mouseButton.button == sf::Mouse::Left)
    {
        sf::FloatRect rect(m_pos.x, m_pos.y, m_size.x, m_size.y);
        if (rect.contains(mp)) m_clicked = true;
    }
}

void Button::draw(sf::RenderTarget& target, const sf::Font& font) const {
    sf::RectangleShape box(m_size);
    box.setPosition(m_pos);
    if (m_active) {
        box.setFillColor(sf::Color(60, 110, 160));
        box.setOutlineColor(sf::Color(150, 200, 255));
    } else {
        box.setFillColor(sf::Color(45, 55, 75));
        box.setOutlineColor(sf::Color(110, 130, 160));
    }
    box.setOutlineThickness(1.5f);
    target.draw(box);

    sf::Text t = makeText(m_label, font, 12, sf::Color(230, 235, 245), {0.f, 0.f});
    const auto bb = t.getLocalBounds();
    t.setPosition(m_pos.x + (m_size.x - bb.width)  * 0.5f - bb.left,
                  m_pos.y + (m_size.y - bb.height) * 0.5f - bb.top);
    target.draw(t);
}


// =============================================================================
// UIPanel layout  (panel size expected: 1360 x 260)
// =============================================================================
UIPanel::UIPanel(const sf::Vector2f& pos, const sf::Vector2f& size)
    : m_pos(pos), m_size(size)
{
    // -------------------- Sliders, row A (top) -----------------------------
    m_tempSlider = std::make_unique<Slider>(
        sf::Vector2f{pos.x + 140.f, pos.y + 60.f},
        240.f, 100.f, 600.f, 300.f, /*log=*/false);

    m_lambdaSlider = std::make_unique<Slider>(
        sf::Vector2f{pos.x + 500.f, pos.y + 60.f},
        240.f, 300.f, 1500.f, 1500.f, /*log=*/false);

    // -------------------- Sliders, row B (below) ---------------------------
    m_dopingSlider = std::make_unique<Slider>(
        sf::Vector2f{pos.x + 140.f, pos.y + 125.f},
        240.f, 1e14f, 1e19f, 1e16f, /*log=*/true);

    m_bFieldSlider = std::make_unique<Slider>(
        sf::Vector2f{pos.x + 500.f, pos.y + 125.f},
        240.f, 0.f, 10.f, 0.f, /*log=*/false);

    // -------------------- Buttons, rows at y=+30 and y=+76 -----------------
    const float bx = pos.x + 800.f;

    m_btnIntrinsic = std::make_unique<Button>(
        sf::Vector2f{bx,         pos.y + 30.f},
        sf::Vector2f{110.f, 32.f}, "Intrinsic");
    m_btnNType = std::make_unique<Button>(
        sf::Vector2f{bx + 120.f, pos.y + 30.f},
        sf::Vector2f{110.f, 32.f}, "n-type (P)");
    m_btnPType = std::make_unique<Button>(
        sf::Vector2f{bx + 240.f, pos.y + 30.f},
        sf::Vector2f{110.f, 32.f}, "p-type (B)");
    m_btnIntrinsic->setActive(true);

    m_btnIonization = std::make_unique<Button>(
        sf::Vector2f{bx,         pos.y + 72.f},
        sf::Vector2f{170.f, 32.f}, "Ionization: Full");
    m_btnLight = std::make_unique<Button>(
        sf::Vector2f{bx + 180.f, pos.y + 72.f},
        sf::Vector2f{170.f, 32.f}, "Light: Off");

    m_btnMobility = std::make_unique<Button>(
        sf::Vector2f{bx,         pos.y + 114.f},
        sf::Vector2f{220.f, 32.f}, "Mobility: Matthiessen");
    m_btnMobility->setActive(true);
    m_btnExportCSV = std::make_unique<Button>(
        sf::Vector2f{bx + 230.f, pos.y + 114.f},
        sf::Vector2f{120.f, 32.f}, "Export CSV");
}


// =============================================================================
// Event dispatch
// =============================================================================
bool UIPanel::handleEvent(const sf::Event&    ev,
                          const sf::Vector2f& mp,
                          PhysicsEngine&      physics)
{
    bool rebuild = false;

    // ---- Widgets ----------------------------------------------------------
    m_tempSlider   ->handleEvent(ev, mp);
    m_dopingSlider ->handleEvent(ev, mp);
    m_lambdaSlider ->handleEvent(ev, mp);
    m_bFieldSlider ->handleEvent(ev, mp);

    m_btnIntrinsic ->handleEvent(ev, mp);
    m_btnNType     ->handleEvent(ev, mp);
    m_btnPType     ->handleEvent(ev, mp);
    m_btnIonization->handleEvent(ev, mp);
    m_btnLight     ->handleEvent(ev, mp);
    m_btnMobility  ->handleEvent(ev, mp);
    m_btnExportCSV ->handleEvent(ev, mp);

    // ---- Sliders ----------------------------------------------------------
    if (m_tempSlider->consumeChanged()) {
        physics.setTemperature(m_tempSlider->getValue());
        rebuild = true;
    }
    if (m_dopingSlider->consumeChanged()) {
        physics.setDopingConcentration(m_dopingSlider->getValue());
        rebuild = true;
    }
    if (m_lambdaSlider->consumeChanged()) {
        physics.setWavelengthNm(m_lambdaSlider->getValue());
        rebuild = true;   // Delta n may change -> visible carrier cloud
    }
    if (m_bFieldSlider->consumeChanged()) {
        physics.setMagneticField(m_bFieldSlider->getValue());
        // No rebuild needed -- carriers keep their state, only dynamics change.
    }

    // ---- Doping type ------------------------------------------------------
    if (m_btnIntrinsic->consumeClicked()) {
        physics.setDopingType(DopingType::Intrinsic);
        m_btnIntrinsic->setActive(true);
        m_btnNType    ->setActive(false);
        m_btnPType    ->setActive(false);
        rebuild = true;
    }
    if (m_btnNType->consumeClicked()) {
        physics.setDopingType(DopingType::NType);
        m_btnIntrinsic->setActive(false);
        m_btnNType    ->setActive(true);
        m_btnPType    ->setActive(false);
        rebuild = true;
    }
    if (m_btnPType->consumeClicked()) {
        physics.setDopingType(DopingType::PType);
        m_btnIntrinsic->setActive(false);
        m_btnNType    ->setActive(false);
        m_btnPType    ->setActive(true);
        rebuild = true;
    }

    // ---- Ionization mode --------------------------------------------------
    if (m_btnIonization->consumeClicked()) {
        const bool newState = !physics.getIncompleteIonization();
        physics.setIncompleteIonization(newState);
        m_btnIonization->setActive(newState);
        m_btnIonization->setLabel(newState ? "Ionization: Incomplete"
                                           : "Ionization: Full");
        rebuild = true;
    }

    // ---- Light toggle -----------------------------------------------------
    if (m_btnLight->consumeClicked()) {
        const bool newState = !physics.getOpticalEnabled();
        physics.setOpticalEnabled(newState);
        m_btnLight->setActive(newState);
        m_btnLight->setLabel(newState ? "Light: On" : "Light: Off");
        rebuild = true;
    }

    // ---- Mobility model ---------------------------------------------------
    if (m_btnMobility->consumeClicked()) {
        const bool toMatthiessen =
            (physics.getMobilityModel() != MobilityModel::Matthiessen);
        physics.setMobilityModel(toMatthiessen
            ? MobilityModel::Matthiessen
            : MobilityModel::Arora);
        m_btnMobility->setLabel(toMatthiessen
            ? "Mobility: Matthiessen"
            : "Mobility: Arora");
        m_btnMobility->setActive(true);
    }

    // ---- CSV export -------------------------------------------------------
    if (m_btnExportCSV->consumeClicked()) {
        performCSVExport(physics);
    }

    return rebuild;
}


// =============================================================================
// CSV export
// =============================================================================
void UIPanel::performCSVExport(const PhysicsEngine& physics) {
    const std::string path = "export_data.csv";
    const bool ok = physics.exportCSV(path);
    m_statusMessage = ok ? "Saved row to " + path
                         : "Failed to write " + path;
    m_statusLifetime = 3.0f;
    m_statusClock.restart();
}


// =============================================================================
// Rendering
// =============================================================================
void UIPanel::draw(sf::RenderTarget&   target,
                   const PhysicsEngine& physics,
                   const sf::Font&      font) const
{
    // Panel background
    sf::RectangleShape bg(m_size);
    bg.setPosition(m_pos);
    bg.setFillColor(sf::Color(28, 33, 45));
    bg.setOutlineColor(sf::Color(80, 90, 110));
    bg.setOutlineThickness(1.f);
    target.draw(bg);

    target.draw(makeText("Controls", font, 16, sf::Color(230, 235, 245),
                         {m_pos.x + 14.f, m_pos.y + 8.f}));

    // --- Slider values strings --------------------------------------------
    std::ostringstream sT, sN, sL, sB;
    sT << std::fixed << std::setprecision(1) << m_tempSlider->getValue() << " K";
    sN << formatSci(m_dopingSlider->getValue()) << " cm^-3";
    sL << std::fixed << std::setprecision(0) << m_lambdaSlider->getValue() << " nm";
    sB << std::fixed << std::setprecision(2) << m_bFieldSlider->getValue() << " T";

    m_tempSlider  ->draw(target, "Temperature (T)",          sT.str(), font);
    m_dopingSlider->draw(target, "Doping (N_d / N_a)",       sN.str(), font);
    m_lambdaSlider->draw(target, "Wavelength (lambda)",      sL.str(), font);
    m_bFieldSlider->draw(target, "Magnetic field (B)",       sB.str(), font);

    // --- Wavelength colour indicator (little swatch next to the slider) ---
    {
        sf::RectangleShape swatch({18.f, 18.f});
        swatch.setPosition(m_pos.x + 500.f + 240.f + 70.f, m_pos.y + 52.f);
        swatch.setFillColor(wavelengthToColor(m_lambdaSlider->getValue()));
        swatch.setOutlineColor(sf::Color(200, 210, 230));
        swatch.setOutlineThickness(1.f);
        target.draw(swatch);
    }

    // --- Buttons ----------------------------------------------------------
    m_btnIntrinsic ->draw(target, font);
    m_btnNType     ->draw(target, font);
    m_btnPType     ->draw(target, font);
    m_btnIonization->draw(target, font);
    m_btnLight     ->draw(target, font);
    m_btnMobility  ->draw(target, font);
    m_btnExportCSV ->draw(target, font);

    // --- Status message ---------------------------------------------------
    if (!m_statusMessage.empty() && m_statusLifetime > 0.f) {
        const float elapsed   = m_statusClock.getElapsedTime().asSeconds();
        const float remaining = m_statusLifetime - elapsed;
        if (remaining > 0.f) {
            const float alpha = std::clamp(remaining / m_statusLifetime, 0.f, 1.f);
            const sf::Color col(180, 220, 255,
                                static_cast<sf::Uint8>(255 * alpha));
            target.draw(makeText(m_statusMessage, font, 12, col,
                                 {m_pos.x + 800.f, m_pos.y + 158.f}));
        }
    }

    // =======================================================================
    // Readouts column (right-hand side)
    // =======================================================================
    const float rX  = m_pos.x + 1180.f;
    float       rY  = m_pos.y + 14.f;
    const auto  lbl = sf::Color(200, 205, 220);

    target.draw(makeText("Physics Readouts", font, 13,
                         sf::Color(230, 235, 245), {rX, rY}));
    rY += 18.f;

    auto row = [&](const std::string& label,
                   const std::string& value,
                   sf::Color valueCol = sf::Color(255, 240, 170))
    {
        target.draw(makeText(label, font, 11, lbl,      {rX,         rY}));
        target.draw(makeText(value, font, 11, valueCol, {rX + 70.f,  rY}));
        rY += 14.f;
    };

    row("T",       formatFixed(physics.getTemperature(), 1) + " K");
    row("E_g",     formatFixed(physics.getBandgap(),   4) + " eV");
    row("E_f",     formatFixed(physics.getFermiLevel(),4) + " eV");
    row("E_phot",  formatFixed(physics.getPhotonEnergy(), 3) + " eV");
    row("n_i",     formatSci  (physics.getIntrinsicCarrier()) + " cm-3");
    row("n",       formatSci  (physics.getTotalElectronConc()) + " cm-3",
        sf::Color(255, 230, 80));
    row("p",       formatSci  (physics.getTotalHoleConc())     + " cm-3",
        sf::Color(255, 100, 200));
    row("dN_opt",  formatSci  (physics.getExcessCarrierDensity()) + " cm-3",
        sf::Color(140, 230, 255));
    row("mu_n",    formatFixed(physics.getElectronMobility(), 1) + " cm2/Vs");
    row("mu_p",    formatFixed(physics.getHoleMobility(),     1) + " cm2/Vs");
    row("sigma",   formatSci  (physics.getConductivity())   + " S/cm",
        sf::Color(140, 255, 180));
    row("rho",     formatSci  (physics.getResistivity())    + " Ohm.cm");
    row("R_H",     formatSci  (physics.getHallCoefficient()) + " cm3/C",
        sf::Color(255, 180, 120));
    row("B",       formatFixed(physics.getMagneticField(), 2) + " T");
    row("Ion %",   formatFixed(physics.getIonizationFraction() * 100.0, 1) + " %");
}
