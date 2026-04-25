#include "UIPanel.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <sstream>

#include "Palette.hpp"

using palette::makeText;
using palette::wavelengthColor;


// =============================================================================
// Local numeric formatting
// =============================================================================
namespace {

    [[nodiscard]] std::string sci(double v, int prec = 2) {
        std::ostringstream oss;
        oss << std::scientific << std::setprecision(prec) << v;
        return oss.str();
    }

    [[nodiscard]] std::string fix(double v, int prec = 2) {
        std::ostringstream oss;
        oss << std::fixed << std::setprecision(prec) << v;
        return oss.str();
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
    return sf::FloatRect({m_pos.x - 6.f, m_pos.y - 12.f},
                         {m_length + 12.f, 24.f});
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
    if (const auto* e = ev.getIf<sf::Event::MouseButtonPressed>()) {
        if (e->button == sf::Mouse::Button::Left && bounds().contains(mp)) {
            m_dragging = true;
            setFromNormalized((mp.x - m_pos.x) / m_length);
        }
    }
    else if (const auto* e = ev.getIf<sf::Event::MouseButtonReleased>()) {
        if (e->button == sf::Mouse::Button::Left) m_dragging = false;
    }
    else if (ev.is<sf::Event::MouseMoved>() && m_dragging) {
        setFromNormalized((mp.x - m_pos.x) / m_length);
    }
}

void Slider::draw(sf::RenderTarget&  target,
                  const std::string& label,
                  const std::string& valueStr,
                  const sf::Font&    font) const
{
    sf::RectangleShape track({m_length, 4.f});
    track.setPosition({m_pos.x, m_pos.y - 2.f});
    track.setFillColor(sf::Color(90, 100, 120));
    target.draw(track);

    const float handleX = m_pos.x + normalized() * m_length;
    sf::CircleShape handle(8.f);
    handle.setOrigin({8.f, 8.f});
    handle.setPosition({handleX, m_pos.y});
    handle.setFillColor(palette::FermiLine);
    handle.setOutlineColor(sf::Color(255, 255, 200));
    handle.setOutlineThickness(1.5f);
    target.draw(handle);

    target.draw(makeText(font, label, 13, palette::TextDim,
                         {m_pos.x, m_pos.y - 28.f}));
    target.draw(makeText(font, valueStr, 13, palette::Accent,
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
    if (const auto* e = ev.getIf<sf::Event::MouseButtonPressed>()) {
        if (e->button == sf::Mouse::Button::Left) {
            const sf::FloatRect rect(m_pos, m_size);
            if (rect.contains(mp)) m_clicked = true;
        }
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

    sf::Text t = makeText(font, m_label, 12, palette::TextLight, {0.f, 0.f});
    const auto bb = t.getLocalBounds();
    t.setPosition({m_pos.x + (m_size.x - bb.size.x) * 0.5f - bb.position.x,
                   m_pos.y + (m_size.y - bb.size.y) * 0.5f - bb.position.y});
    target.draw(t);
}


// =============================================================================
// UIPanel layout (panel size: 1360 x 260)
// =============================================================================
UIPanel::UIPanel(const sf::Vector2f& pos, const sf::Vector2f& size)
    : m_pos(pos), m_size(size)
{
    // -------------------- Sliders -----------------------------------------
    m_tempSlider = std::make_unique<Slider>(
        sf::Vector2f{pos.x + 140.f, pos.y + 60.f},
        240.f, 100.f, 600.f, 300.f, /*log=*/false);

    m_lambdaSlider = std::make_unique<Slider>(
        sf::Vector2f{pos.x + 500.f, pos.y + 60.f},
        240.f, 300.f, 1500.f, 1500.f, /*log=*/false);

    m_dopingSlider = std::make_unique<Slider>(
        sf::Vector2f{pos.x + 140.f, pos.y + 125.f},
        240.f, 1e14f, 1e19f, 1e16f, /*log=*/true);

    m_bFieldSlider = std::make_unique<Slider>(
        sf::Vector2f{pos.x + 500.f, pos.y + 125.f},
        240.f, 0.f, 10.f, 0.f, /*log=*/false);

    // -------------------- Buttons -----------------------------------------
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

    // Forward event to every widget.
    const std::array<Slider*, 4> sliders{
        m_tempSlider.get(), m_dopingSlider.get(),
        m_lambdaSlider.get(), m_bFieldSlider.get()
    };
    for (auto* s : sliders) s->handleEvent(ev, mp);

    const std::array<Button*, 7> buttons{
        m_btnIntrinsic.get(), m_btnNType.get(), m_btnPType.get(),
        m_btnIonization.get(), m_btnLight.get(), m_btnMobility.get(),
        m_btnExportCSV.get()
    };
    for (auto* b : buttons) b->handleEvent(ev, mp);

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
        rebuild = true;
    }
    if (m_bFieldSlider->consumeChanged()) {
        physics.setMagneticField(m_bFieldSlider->getValue());
    }

    // ---- Doping-type radio ------------------------------------------------
    auto selectDopingType = [&](DopingType dt, Button* pressed) {
        physics.setDopingType(dt);
        m_btnIntrinsic->setActive(pressed == m_btnIntrinsic.get());
        m_btnNType    ->setActive(pressed == m_btnNType.get());
        m_btnPType    ->setActive(pressed == m_btnPType.get());
        rebuild = true;
    };
    if (m_btnIntrinsic->consumeClicked())
        selectDopingType(DopingType::Intrinsic, m_btnIntrinsic.get());
    if (m_btnNType->consumeClicked())
        selectDopingType(DopingType::NType, m_btnNType.get());
    if (m_btnPType->consumeClicked())
        selectDopingType(DopingType::PType, m_btnPType.get());

    // ---- Toggles ----------------------------------------------------------
    if (m_btnIonization->consumeClicked()) {
        const bool on = !physics.getIncompleteIonization();
        physics.setIncompleteIonization(on);
        m_btnIonization->setActive(on);
        m_btnIonization->setLabel(on ? "Ionization: Incomplete"
                                     : "Ionization: Full");
        rebuild = true;
    }

    if (m_btnLight->consumeClicked()) {
        const bool on = !physics.getOpticalEnabled();
        physics.setOpticalEnabled(on);
        m_btnLight->setActive(on);
        m_btnLight->setLabel(on ? "Light: On" : "Light: Off");
        rebuild = true;
    }

    if (m_btnMobility->consumeClicked()) {
        const bool toMatthiessen =
            (physics.getMobilityModel() != MobilityModel::Matthiessen);
        physics.setMobilityModel(toMatthiessen
            ? MobilityModel::Matthiessen
            : MobilityModel::Arora);
        m_btnMobility->setLabel(toMatthiessen ? "Mobility: Matthiessen"
                                              : "Mobility: Arora");
    }

    // ---- Action -----------------------------------------------------------
    if (m_btnExportCSV->consumeClicked()) performCSVExport(physics);

    return rebuild;
}


// =============================================================================
// CSV export
// =============================================================================
void UIPanel::performCSVExport(const PhysicsEngine& physics) {
    constexpr const char* path = "export_data.csv";
    const bool ok = physics.exportCSV(path);
    m_statusMessage  = ok ? std::string{"Saved row to "}   + path
                          : std::string{"Failed to write "} + path;
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
    bg.setFillColor(palette::PanelBg);
    bg.setOutlineColor(palette::PanelEdge);
    bg.setOutlineThickness(1.f);
    target.draw(bg);

    target.draw(makeText(font, "Controls", 16, palette::TextLight,
                         {m_pos.x + 14.f, m_pos.y + 8.f}));

    // Slider value strings
    m_tempSlider->draw(target,  "Temperature (T)",
        fix(m_tempSlider->getValue(), 1) + " K", font);
    m_dopingSlider->draw(target, "Doping (N_d / N_a)",
        sci(m_dopingSlider->getValue()) + " cm^-3", font);
    m_lambdaSlider->draw(target, "Wavelength (lambda)",
        fix(m_lambdaSlider->getValue(), 0) + " nm", font);
    m_bFieldSlider->draw(target, "Magnetic field (B)",
        fix(m_bFieldSlider->getValue(), 2) + " T", font);

    // Wavelength colour swatch
    {
        sf::RectangleShape swatch({18.f, 18.f});
        swatch.setPosition({m_pos.x + 500.f + 240.f + 70.f, m_pos.y + 52.f});
        swatch.setFillColor(wavelengthColor(m_lambdaSlider->getValue()));
        swatch.setOutlineColor(sf::Color(200, 210, 230));
        swatch.setOutlineThickness(1.f);
        target.draw(swatch);
    }

    // Buttons
    m_btnIntrinsic ->draw(target, font);
    m_btnNType     ->draw(target, font);
    m_btnPType     ->draw(target, font);
    m_btnIonization->draw(target, font);
    m_btnLight     ->draw(target, font);
    m_btnMobility  ->draw(target, font);
    m_btnExportCSV ->draw(target, font);

    // Fading status message
    if (!m_statusMessage.empty() && m_statusLifetime > 0.f) {
        const float elapsed   = m_statusClock.getElapsedTime().asSeconds();
        const float remaining = m_statusLifetime - elapsed;
        if (remaining > 0.f) {
            const float alpha = std::clamp(remaining / m_statusLifetime, 0.f, 1.f);
            const sf::Color col(180, 220, 255,
                                static_cast<std::uint8_t>(255 * alpha));
            target.draw(makeText(font, m_statusMessage, 12, col,
                                 {m_pos.x + 800.f, m_pos.y + 158.f}));
        }
    }

    drawReadouts(target, physics, font);
}


// =============================================================================
// Readouts column (right-hand side of the panel)
// =============================================================================
void UIPanel::drawReadouts(sf::RenderTarget&   target,
                           const PhysicsEngine& physics,
                           const sf::Font&      font) const
{
    const float rX = m_pos.x + 1170.f;
    float       rY = m_pos.y + 14.f;

    target.draw(makeText(font, "Physics Readouts", 13,
                         palette::TextLight, {rX, rY}));
    rY += 18.f;

    auto row = [&](const std::string& label,
                   const std::string& value,
                   sf::Color           valueCol = palette::Accent)
    {
        target.draw(makeText(font, label, 11, palette::TextDim,
                             {rX, rY}));
        target.draw(makeText(font, value, 11, valueCol,
                             {rX + 80.f, rY}));
        rY += 14.f;
    };

    // --- Thermodynamics ---
    row("T",      fix(physics.getTemperature(), 1) + " K");
    row("E_g",    fix(physics.getBandgap(),     4) + " eV");
    row("E_f",    fix(physics.getFermiLevel(),  4) + " eV");
    row("E_phot", fix(physics.getPhotonEnergy(),3) + " eV");
    row("n_i",    sci(physics.getIntrinsicCarrier())    + " cm-3");

    // --- Carriers ---
    row("n",      sci(physics.getTotalElectronConc())   + " cm-3",
        palette::Electron);
    row("p",      sci(physics.getTotalHoleConc())       + " cm-3",
        palette::Hole);
    if (physics.isOpticallyPumped()) {
        row("dN",  sci(physics.getExcessCarrierDensity()) + " cm-3",
            palette::Optical);
    }

    // --- Transport ---
    row("mu_n",   fix(physics.getElectronMobility(), 1) + " cm2/Vs");
    row("mu_p",   fix(physics.getHoleMobility(),     1) + " cm2/Vs");
    row("sigma",  sci(physics.getConductivity())         + " S/cm",
        palette::Conductivity);
    row("rho",    sci(physics.getResistivity())          + " Ohm.cm");

    // --- Magnetic / Hall ---
    row("R_H",    sci(physics.getHallCoefficient()) + " cm3/C",
        palette::Hall);

    // --- Ionization (only meaningful when doped) ---
    if (physics.getDopingType() != DopingType::Intrinsic) {
        row("Ion%", fix(physics.getIonizationFraction() * 100.0, 1) + " %");
    }
}
