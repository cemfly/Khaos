// =============================================================================
// main.cpp -- Semiconductor Analysis & Simulation Platform
//
//   Author : dex / cemfly-april2026
//   License: MIT
// -----------------------------------------------------------------------------
// Window, ImGui dock-space, panel composition, event routing.
//
// Owned components:
//
//   PhysicsEngine     -- material-aware semiconductor physics core
//   DriftDiffusion    -- 2D FTCS solver for click-deposited sources
//   CrystalView       -- RenderTexture: lattice, carriers, heatmap, B vectors
//   BandView          -- RenderTexture: band diagram + Fermi-Dirac curve
//
// Six dockable ImGui windows:
//   * Controls               -- material picker, sliders, mode toggles
//   * Readouts               -- live numerical state
//   * Live Oscilloscope      -- ImPlot: sigma(t), n(t), p(t)
//   * Crystal View           -- ImGui::Image of CrystalView's RenderTexture
//   * Band Diagram           -- ImGui::Image of BandView's RenderTexture
//   * Drift-Diffusion        -- click controls, source intensity, clear
//
// =============================================================================

#include <SFML/Graphics.hpp>

#include <imgui.h>
#include <imgui-SFML.h>
#include <imgui_internal.h>      // DockBuilder API for the fixed layout
#include <implot.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <filesystem>
#include <iostream>
#include <memory>
#include <optional>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include "BandView.hpp"
#include "CrystalView.hpp"
#include "DriftDiffusion.hpp"
#include "Material.hpp"
#include "Palette.hpp"
#include "PhysicsEngine.hpp"


// =============================================================================
// Configuration
// =============================================================================
namespace {

constexpr unsigned int kWindowWidth  = 1700;
constexpr unsigned int kWindowHeight = 1000;
constexpr const char*  kWindowTitle  =
    "Semiconductor Analysis & Simulation Platform";

constexpr std::size_t  kOscilloscopeMaxSamples = 600;


// -----------------------------------------------------------------------------
// Font loading (filesystem-aware so we don't spam SFML warnings)
// -----------------------------------------------------------------------------
constexpr std::array<std::string_view, 7> kFontCandidates{
    "assets/font.ttf",
    "../assets/font.ttf",
    "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf",
    "/usr/share/fonts/truetype/liberation/LiberationSans-Regular.ttf",
    "/usr/share/fonts/TTF/DejaVuSans.ttf",
    "C:/Windows/Fonts/arial.ttf",
    "C:/Windows/Fonts/segoeui.ttf",
};

[[nodiscard]] bool loadFont(sf::Font& font) {
    for (const auto path : kFontCandidates) {
        const std::filesystem::path fp{path};
        std::error_code ec;
        if (!std::filesystem::exists(fp, ec)) continue;
        if (font.openFromFile(fp)) {
            std::cout << "[info] Loaded font: " << fp.string() << '\n';
            return true;
        }
    }
    std::cerr << "[warn] No usable font found.\n";
    return false;
}


// =============================================================================
// ImGui setup helpers
// =============================================================================
void enableImGuiDocking() {
    ImGuiIO& io = ImGui::GetIO();
    io.ConfigFlags |= ImGuiConfigFlags_DockingEnable;
    io.ConfigFlags |= ImGuiConfigFlags_NavEnableKeyboard;
}

// =============================================================================
// applyLightStyle -- light "scientific" theme, COMSOL / MATLAB / Mathematica
// inspired. White-ish backgrounds, dark slate text, mid-blue accents.
// =============================================================================
void applyLightStyle() {
    ImGui::StyleColorsLight();
    ImGuiStyle& s = ImGui::GetStyle();

    // Roundings -- subtle, no Material-Design exaggeration.
    s.WindowRounding    = 5.0f;
    s.ChildRounding     = 4.0f;
    s.FrameRounding     = 3.0f;
    s.GrabRounding      = 3.0f;
    s.PopupRounding     = 4.0f;
    s.TabRounding       = 3.0f;

    s.WindowBorderSize  = 1.0f;
    s.ChildBorderSize   = 1.0f;
    s.FrameBorderSize   = 1.0f;
    s.PopupBorderSize   = 1.0f;

    // Generous spacing -- breathing room.
    s.WindowPadding     = ImVec2(12.0f, 12.0f);
    s.FramePadding      = ImVec2( 8.0f,  5.0f);
    s.CellPadding       = ImVec2( 6.0f,  4.0f);
    s.ItemSpacing       = ImVec2(10.0f,  8.0f);
    s.ItemInnerSpacing  = ImVec2( 6.0f,  6.0f);
    s.IndentSpacing     = 22.0f;
    s.ScrollbarSize     = 14.0f;
    s.GrabMinSize       = 12.0f;

    ImVec4* c = s.Colors;

    // Backgrounds  (white  ->  near-white  ->  off-white)
    c[ImGuiCol_WindowBg]            = ImVec4(0.972f, 0.976f, 0.984f, 1.00f);
    c[ImGuiCol_ChildBg]             = ImVec4(1.000f, 1.000f, 1.000f, 1.00f);
    c[ImGuiCol_PopupBg]             = ImVec4(0.985f, 0.985f, 0.990f, 0.98f);
    c[ImGuiCol_MenuBarBg]           = ImVec4(0.945f, 0.950f, 0.960f, 1.00f);

    // Borders / separators (cool grey)
    c[ImGuiCol_Border]              = ImVec4(0.725f, 0.764f, 0.823f, 0.65f);
    c[ImGuiCol_BorderShadow]        = ImVec4(0.000f, 0.000f, 0.000f, 0.00f);
    c[ImGuiCol_Separator]           = ImVec4(0.745f, 0.784f, 0.843f, 0.65f);
    c[ImGuiCol_SeparatorHovered]    = ImVec4(0.235f, 0.431f, 0.784f, 0.78f);
    c[ImGuiCol_SeparatorActive]     = ImVec4(0.235f, 0.431f, 0.784f, 1.00f);

    // Text
    c[ImGuiCol_Text]                = ImVec4(0.110f, 0.133f, 0.180f, 1.00f);
    c[ImGuiCol_TextDisabled]        = ImVec4(0.470f, 0.510f, 0.580f, 1.00f);

    // Frame (input bg)
    c[ImGuiCol_FrameBg]             = ImVec4(0.960f, 0.965f, 0.975f, 1.00f);
    c[ImGuiCol_FrameBgHovered]      = ImVec4(0.890f, 0.910f, 0.945f, 1.00f);
    c[ImGuiCol_FrameBgActive]       = ImVec4(0.815f, 0.870f, 0.945f, 1.00f);

    // Title bars
    c[ImGuiCol_TitleBg]             = ImVec4(0.910f, 0.925f, 0.945f, 1.00f);
    c[ImGuiCol_TitleBgActive]       = ImVec4(0.745f, 0.835f, 0.945f, 1.00f);
    c[ImGuiCol_TitleBgCollapsed]    = ImVec4(0.910f, 0.925f, 0.945f, 0.75f);

    // CollapsingHeader
    c[ImGuiCol_Header]              = ImVec4(0.745f, 0.835f, 0.945f, 0.65f);
    c[ImGuiCol_HeaderHovered]       = ImVec4(0.470f, 0.670f, 0.910f, 0.78f);
    c[ImGuiCol_HeaderActive]        = ImVec4(0.235f, 0.510f, 0.870f, 0.85f);

    // Buttons
    c[ImGuiCol_Button]              = ImVec4(0.870f, 0.910f, 0.965f, 1.00f);
    c[ImGuiCol_ButtonHovered]       = ImVec4(0.700f, 0.815f, 0.945f, 1.00f);
    c[ImGuiCol_ButtonActive]        = ImVec4(0.470f, 0.670f, 0.910f, 1.00f);

    // Sliders
    c[ImGuiCol_SliderGrab]          = ImVec4(0.235f, 0.431f, 0.784f, 1.00f);
    c[ImGuiCol_SliderGrabActive]    = ImVec4(0.137f, 0.333f, 0.686f, 1.00f);

    // Tabs
    c[ImGuiCol_Tab]                 = ImVec4(0.890f, 0.910f, 0.945f, 1.00f);
    c[ImGuiCol_TabHovered]          = ImVec4(0.700f, 0.815f, 0.945f, 1.00f);
    c[ImGuiCol_TabActive]           = ImVec4(0.470f, 0.670f, 0.910f, 1.00f);
    c[ImGuiCol_TabUnfocused]        = ImVec4(0.910f, 0.925f, 0.945f, 1.00f);
    c[ImGuiCol_TabUnfocusedActive]  = ImVec4(0.815f, 0.870f, 0.945f, 1.00f);

    // Docking
    c[ImGuiCol_DockingPreview]      = ImVec4(0.235f, 0.510f, 0.870f, 0.55f);
    c[ImGuiCol_DockingEmptyBg]      = ImVec4(0.945f, 0.950f, 0.960f, 1.00f);

    // Plots
    c[ImGuiCol_PlotLines]           = ImVec4(0.235f, 0.431f, 0.784f, 1.00f);
    c[ImGuiCol_PlotLinesHovered]    = ImVec4(0.870f, 0.215f, 0.215f, 1.00f);
    c[ImGuiCol_PlotHistogram]       = ImVec4(0.235f, 0.431f, 0.784f, 1.00f);
}


// =============================================================================
// HelpMarker -- the standard ImGui "(?)" hover tooltip pattern.
//
// Replaces the wall-of-text TextWrapped() blocks that used to sit inside
// each panel. Now every formula / explanation is opt-in: the user only sees
// it if they hover the icon.
// =============================================================================
void HelpMarker(const char* desc) {
    ImGui::SameLine();
    ImGui::TextDisabled("(?)");
    if (ImGui::BeginItemTooltip()) {
        ImGui::PushTextWrapPos(ImGui::GetFontSize() * 32.0f);
        ImGui::TextUnformatted(desc);
        ImGui::PopTextWrapPos();
        ImGui::EndTooltip();
    }
}


// =============================================================================
// Initial dock layout -- 3-column IDE arrangement.
//
//   Left  25%  : Controls (single tall panel)
//   Centre 50% : top -> Crystal View, bottom -> [Live Osc | Spectrum | I-V] tabs
//   Right 25%  : top -> Band Diagram,  bottom -> [Readouts | Crystal Info] tabs
//
// Called only when no layout exists yet (first launch) or when the user picks
// "View -> Reset layout".
// =============================================================================
void setupInitialDockLayout(ImGuiID dockspace_id, ImVec2 size) {
    ImGui::DockBuilderRemoveNode(dockspace_id);
    ImGui::DockBuilderAddNode(dockspace_id, ImGuiDockNodeFlags_DockSpace);
    ImGui::DockBuilderSetNodeSize(dockspace_id, size);

    ImGuiID dock_main = dockspace_id;

    // 1) Carve off the left 25% for Controls.
    const ImGuiID dock_left =
        ImGui::DockBuilderSplitNode(dock_main, ImGuiDir_Left,
                                    0.25f, nullptr, &dock_main);

    // 2) Carve off the right column (1/3 of the remaining 75% = 25% total).
    ImGuiID dock_right =
        ImGui::DockBuilderSplitNode(dock_main, ImGuiDir_Right,
                                    0.333f, nullptr, &dock_main);

    // 3) Split the centre column 60% top / 40% bottom.
    const ImGuiID dock_center_bottom =
        ImGui::DockBuilderSplitNode(dock_main, ImGuiDir_Down,
                                    0.40f, nullptr, &dock_main);
    const ImGuiID dock_center_top = dock_main;

    // 4) Split the right column 50% top / 50% bottom.
    const ImGuiID dock_right_bottom =
        ImGui::DockBuilderSplitNode(dock_right, ImGuiDir_Down,
                                    0.50f, nullptr, &dock_right);
    const ImGuiID dock_right_top = dock_right;

    // 5) Dock every named window into its slot.  Multiple windows in the
    //    same slot become tabs automatically.
    ImGui::DockBuilderDockWindow("Controls",          dock_left);
    ImGui::DockBuilderDockWindow("Crystal View",      dock_center_top);
    ImGui::DockBuilderDockWindow("Live Oscilloscope", dock_center_bottom);
    ImGui::DockBuilderDockWindow("Spectrum",          dock_center_bottom);
    ImGui::DockBuilderDockWindow("I-V Curve",         dock_center_bottom);
    ImGui::DockBuilderDockWindow("Band Diagram",      dock_right_top);
    ImGui::DockBuilderDockWindow("Readouts",          dock_right_bottom);
    ImGui::DockBuilderDockWindow("Crystal Info",      dock_right_bottom);

    ImGui::DockBuilderFinish(dockspace_id);
}

// Returns the dockspace id so the caller (menu code) can request a layout
// rebuild. The first frame builds the default 3-column layout if no saved
// layout exists.
ImGuiID beginDockspaceHost(bool& requestResetLayout) {
    constexpr ImGuiWindowFlags hostFlags =
          ImGuiWindowFlags_NoTitleBar
        | ImGuiWindowFlags_NoCollapse
        | ImGuiWindowFlags_NoResize
        | ImGuiWindowFlags_NoMove
        | ImGuiWindowFlags_NoBringToFrontOnFocus
        | ImGuiWindowFlags_NoNavFocus
        | ImGuiWindowFlags_MenuBar
        | ImGuiWindowFlags_NoBackground;

    const ImGuiViewport* viewport = ImGui::GetMainViewport();
    ImGui::SetNextWindowPos(viewport->WorkPos);
    ImGui::SetNextWindowSize(viewport->WorkSize);
    ImGui::SetNextWindowViewport(viewport->ID);

    ImGui::PushStyleVar(ImGuiStyleVar_WindowRounding,   0.0f);
    ImGui::PushStyleVar(ImGuiStyleVar_WindowBorderSize, 0.0f);
    ImGui::PushStyleVar(ImGuiStyleVar_WindowPadding,    ImVec2(0.0f, 0.0f));
    ImGui::Begin("##DockHost", nullptr, hostFlags);
    ImGui::PopStyleVar(3);

    const ImGuiID dockspace_id = ImGui::GetID("MainDockSpace");

    // First-run layout build, or explicit reset via menu.
    if (requestResetLayout
        || ImGui::DockBuilderGetNode(dockspace_id) == nullptr)
    {
        setupInitialDockLayout(dockspace_id, viewport->WorkSize);
        requestResetLayout = false;
    }

    ImGui::DockSpace(dockspace_id, ImVec2(0.0f, 0.0f),
                     ImGuiDockNodeFlags_None);

    return dockspace_id;
}


// =============================================================================
// UI state shared across windows
// =============================================================================
struct UIState {
    HeatmapMode heatmapMode    = HeatmapMode::Carriers;
    bool        showVectorField = true;
    bool        ddInteractive   = true;     // click adds drift-diffusion source
    float       sourceIntensity = 1.5f;
    float       sourceSigma     = 0.05f;

    // Phase 6 -- BJT controls
    float       V_BE = 0.0f;
    float       V_CE = 0.0f;

    std::string statusMessage;
    sf::Clock   statusClock;
    float       statusLifetime = 0.0f;

    void flashStatus(std::string s, float life = 3.0f) {
        statusMessage  = std::move(s);
        statusLifetime = life;
        statusClock.restart();
    }
};


// =============================================================================
// Menu bar
// =============================================================================
void drawMenuBar(bool& running, UIState& ui, PhysicsEngine& physics,
                 DriftDiffusion& dd, bool& requestResetLayout)
{
    if (ImGui::BeginMenuBar()) {
        if (ImGui::BeginMenu("File")) {
            if (ImGui::MenuItem("Export CSV row")) {
                if (physics.exportCSV("export_data.csv"))
                    ui.flashStatus("Saved row to export_data.csv");
                else
                    ui.flashStatus("Failed to write export_data.csv");
            }
            ImGui::Separator();
            if (ImGui::MenuItem("Exit", "Esc")) running = false;
            ImGui::EndMenu();
        }
        if (ImGui::BeginMenu("View")) {
            if (ImGui::BeginMenu("Heatmap layer")) {
                bool h_none = ui.heatmapMode == HeatmapMode::None;
                bool h_n    = ui.heatmapMode == HeatmapMode::Carriers;
                bool h_T    = ui.heatmapMode == HeatmapMode::Thermal;
                if (ImGui::MenuItem("Off",             nullptr, &h_none))
                    ui.heatmapMode = HeatmapMode::None;
                if (ImGui::MenuItem("Carriers n(x,y)", nullptr, &h_n))
                    ui.heatmapMode = HeatmapMode::Carriers;
                if (ImGui::MenuItem("Thermal T(x,y)",  nullptr, &h_T))
                    ui.heatmapMode = HeatmapMode::Thermal;
                ImGui::EndMenu();
            }
            ImGui::MenuItem("Lorentz vectors",   nullptr, &ui.showVectorField);
            ImGui::MenuItem("Click adds source", nullptr, &ui.ddInteractive);
            ImGui::Separator();
            if (ImGui::MenuItem("Clear sources / reset thermal")) {
                dd.clear();
                physics.setDriftDiffusionExcess(0.0);
                ui.flashStatus("Drift-diffusion + thermal grid cleared");
            }
            if (ImGui::MenuItem("Reset window layout")) {
                requestResetLayout = true;
                ui.flashStatus("Layout restored to default");
            }
            ImGui::EndMenu();
        }
        if (ImGui::BeginMenu("Help")) {
            ImGui::TextDisabled("Phase 6 -- electrothermal solver + NPN BJT");
            ImGui::EndMenu();
        }
        ImGui::EndMenuBar();
    }

    if (!ui.statusMessage.empty()) {
        const float elapsed = ui.statusClock.getElapsedTime().asSeconds();
        if (elapsed < ui.statusLifetime) {
            ImGui::SameLine();
            ImGui::TextColored(ImVec4(0.7f, 0.9f, 1.0f, 1.0f),
                               "%s", ui.statusMessage.c_str());
        }
    }
}


// =============================================================================
// Controls window
// =============================================================================
void drawControlsWindow(PhysicsEngine& physics,
                        DriftDiffusion& dd,
                        UIState& ui)
{
    if (!ImGui::Begin("Controls")) { ImGui::End(); return; }

    // =====================================================================
    // Material & Doping
    // =====================================================================
    if (ImGui::CollapsingHeader("Material & Doping",
                                ImGuiTreeNodeFlags_DefaultOpen))
    {
        int current = static_cast<int>(physics.getMaterialKind());
        if (ImGui::Combo("Material", &current,
                         material::kLabels, material::kCount))
        {
            physics.setMaterial(static_cast<material::Kind>(current));
            dd.configureForMaterial(physics.getMaterial());
        }
        HelpMarker(
            "Switches the physics engine between Si / GaAs / Ge profiles. "
            "Each one carries its own Varshni coefficients, effective DOS, "
            "Matthiessen mobility constants, optical absorption character "
            "(direct vs indirect) and thermal conductivity.");

        const auto& m = physics.getMaterial();
        ImGui::TextDisabled("%s  |  E_g(0)=%.3f eV  |  %s",
            std::string(m.name).c_str(), m.Eg0,
            m.isDirectBandgap ? "DIRECT bandgap" : "indirect bandgap");

        const char* dopingLabels[] = { "Intrinsic", "n-type", "p-type" };
        int dt_idx = static_cast<int>(physics.getDopingType());
        if (ImGui::Combo("Doping type", &dt_idx, dopingLabels,
                         IM_ARRAYSIZE(dopingLabels)))
            physics.setDopingType(static_cast<DopingType>(dt_idx));
        HelpMarker(
            "Intrinsic: pure host, n = p = n_i. n-type adds donor atoms "
            "(P in Si, Si on Ga site in GaAs). p-type adds acceptors "
            "(B in Si, Be in GaAs).");

        float logN = static_cast<float>(std::log10(
            std::max(physics.getDopingConcentration(), 1.0e10)));
        if (ImGui::SliderFloat("log10 N", &logN, 14.0f, 19.0f, "%.2f"))
            physics.setDopingConcentration(std::pow(10.0f, logN));
        HelpMarker(
            "Doping concentration in cm^-3 on a logarithmic scale. "
            "Slider value is log10 of the dopant density.");
        ImGui::SameLine();
        ImGui::TextDisabled("(%.2e cm^-3)", physics.getDopingConcentration());

        bool ii = physics.getIncompleteIonization();
        if (ImGui::Checkbox("Incomplete ionization (freeze-out)", &ii))
            physics.setIncompleteIonization(ii);
        HelpMarker(
            "When enabled, the donor / acceptor ionization is governed by "
            "Fermi-Dirac statistics:\n"
            "  N_d^+ = N_d / (1 + g_D exp((E_f - E_d)/kT))\n"
            "  N_a^- = N_a / (1 + g_A exp((E_a - E_f)/kT))\n"
            "At cryogenic temperatures, this exposes the carrier freeze-out "
            "regime where n << N_d.");
    }

    // =====================================================================
    // Thermal & Transport
    // =====================================================================
    if (ImGui::CollapsingHeader("Thermal & Transport",
                                ImGuiTreeNodeFlags_DefaultOpen))
    {
        float T = static_cast<float>(physics.getTemperature());
        if (ImGui::SliderFloat("Temperature [K]", &T, 100.0f, 600.0f, "%.1f"))
            physics.setTemperature(T);
        HelpMarker(
            "Lattice temperature. Drives the Varshni bandgap E_g(T), the "
            "effective DOS (T^(3/2) scaling) and the lattice-limited "
            "mobility (T^(-3/2)). Also sets the Dirichlet boundary "
            "condition on the 2D thermal grid.");

        const char* mlabs[] = { "Matthiessen", "Arora" };
        int mm = static_cast<int>(physics.getMobilityModel());
        if (ImGui::Combo("Mobility model", &mm, mlabs, IM_ARRAYSIZE(mlabs)))
            physics.setMobilityModel(static_cast<MobilityModel>(mm));
        HelpMarker(
            "Matthiessen: 1/mu = 1/mu_lattice + 1/mu_impurity, with explicit "
            "T^(-3/2) phonon and T^(3/2)/N impurity scattering scaling.\n"
            "Arora: empirical Caughey-Thomas form fitted on Si experimental "
            "data; useful for cross-checking the Matthiessen prediction.");
    }

    // =====================================================================
    // Optical
    // =====================================================================
    if (ImGui::CollapsingHeader("Optical"))
    {
        bool light = physics.getOpticalEnabled();
        if (ImGui::Checkbox("Light source ON", &light))
            physics.setOpticalEnabled(light);
        HelpMarker(
            "Toggles the photon source. When the photon energy hv = hc/lambda "
            "exceeds E_g, electron-hole pairs are generated and Delta n is "
            "added to the steady-state carrier population. GaAs (direct gap) "
            "shows a much steeper response than Si (indirect).");

        float lambda = static_cast<float>(physics.getWavelengthNm());
        if (ImGui::SliderFloat("Wavelength [nm]", &lambda,
                               300.0f, 1500.0f, "%.0f"))
            physics.setWavelengthNm(lambda);
        ImGui::SameLine();
        ImGui::TextDisabled("(E=%.2f eV)",
            PhysicsEngine::photonEnergyEv(physics.getWavelengthNm()));
    }

    // =====================================================================
    // Magnetic
    // =====================================================================
    if (ImGui::CollapsingHeader("Magnetic"))
    {
        float B = static_cast<float>(physics.getMagneticField());
        if (ImGui::SliderFloat("B field [T]", &B, 0.0f, 10.0f, "%.2f"))
            physics.setMagneticField(B);
        HelpMarker(
            "Magnetic field perpendicular to the screen plane. Drives "
            "the Lorentz force F = q(v x B) on free carriers (visible as "
            "opposite curl directions for electrons vs holes) and sets the "
            "sign / magnitude of the Hall coefficient R_H.");
    }

    // =====================================================================
    // Device & BJT
    // =====================================================================
    if (ImGui::CollapsingHeader("Device & BJT"))
    {
        const char* devLabels[] = { "Bulk wafer", "NPN BJT" };
        int dev = static_cast<int>(dd.deviceMode());
        if (ImGui::Combo("Device mode", &dev, devLabels, IM_ARRAYSIZE(devLabels)))
            dd.setDeviceMode(static_cast<DeviceMode>(dev));
        HelpMarker(
            "Bulk: a single homogeneous wafer; click in the Crystal View to "
            "deposit a Gaussian generation source.\n"
            "NPN BJT: the grid is partitioned into Emitter | Base | Collector. "
            "V_BE forward-biases the base-emitter junction (electron "
            "injection); V_CE reverse-biases the collector and sweeps "
            "carriers across the base.");

        if (dd.deviceMode() == DeviceMode::NpnBjt) {
            if (ImGui::SliderFloat("V_BE [V]", &ui.V_BE, 0.0f, 1.0f, "%.3f"))
                dd.setBjtVoltages(ui.V_BE, ui.V_CE);
            HelpMarker(
                "Base-emitter forward bias. Emitter electron concentration "
                "follows  n_E = n_0 * exp(q V_BE / kT)  --  exponential in "
                "V_BE.");

            if (ImGui::SliderFloat("V_CE [V]", &ui.V_CE, 0.0f, 5.0f, "%.2f"))
                dd.setBjtVoltages(ui.V_BE, ui.V_CE);
            HelpMarker(
                "Collector-emitter bias. Drives the collector sweep-out "
                "rate and contributes to Joule heating along the active "
                "path. Together with V_BE it defines the BJT operating "
                "point (cut-off / forward-active / saturation).");
        }
    }

    // =====================================================================
    // Visualization
    // =====================================================================
    if (ImGui::CollapsingHeader("Visualization"))
    {
        const char* hLabels[] = { "Off", "Carriers n(x,y)", "Thermal T(x,y)" };
        int hm = static_cast<int>(ui.heatmapMode);
        if (ImGui::Combo("Heatmap layer", &hm, hLabels, IM_ARRAYSIZE(hLabels)))
            ui.heatmapMode = static_cast<HeatmapMode>(hm);
        HelpMarker(
            "Off: lattice + carriers only.\n"
            "Carriers: viridis overlay of the n(x,y) grid (laser spots, "
            "BJT injection profile).\n"
            "Thermal: blue (300 K) -> yellow / red / white (>900 K) overlay "
            "of T(x,y). Use this to spot thermal runaway in BJT mode.");

        ImGui::Checkbox("Lorentz vectors",  &ui.showVectorField);
        ImGui::Checkbox("Click adds source",&ui.ddInteractive);
        HelpMarker(
            "When enabled, left-clicking inside the Crystal View deposits a "
            "Gaussian generation source onto the carrier grid. Use this to "
            "explore drift-diffusion in bulk mode.");

        ImGui::SliderFloat("Source intensity", &ui.sourceIntensity, 0.1f, 5.0f);
        ImGui::SliderFloat("Source size",      &ui.sourceSigma,     0.02f, 0.20f);

        if (ImGui::Button("Clear sources / reset thermal grid")) {
            dd.clear();
            physics.setDriftDiffusionExcess(0.0);
        }
    }

    ImGui::End();
}


// =============================================================================
// Readouts window
// =============================================================================
void drawReadoutsWindow(const PhysicsEngine& physics,
                        const DriftDiffusion& dd) {
    if (!ImGui::Begin("Readouts")) { ImGui::End(); return; }

    auto row = [](const char* label, double v, const char* fmt,
                  const ImVec4& col = ImVec4(1, 1, 1, 1))
    {
        ImGui::Text("%-10s", label);
        ImGui::SameLine();
        ImGui::TextColored(col, fmt, v);
    };

    ImGui::SeparatorText("Material");
    ImGui::Text("%-10s %s", "Profile",
        std::string(physics.getMaterial().name).c_str());

    ImGui::SeparatorText("Thermodynamic");
    row("T",      physics.getTemperature(),         "%8.2f K");
    row("E_g",    physics.getBandgap(),             "%8.4f eV");
    row("E_f",    physics.getFermiLevel(),          "%8.4f eV");
    row("n_i",    physics.getIntrinsicCarrier(),    "%.3e cm^-3");

    ImGui::SeparatorText("Carriers (total = thermal + optical + drift)");
    row("n",      physics.getTotalElectronConc(),   "%.3e cm^-3",
        ImVec4(1.0f, 0.9f, 0.3f, 1.0f));
    row("p",      physics.getTotalHoleConc(),       "%.3e cm^-3",
        ImVec4(1.0f, 0.4f, 0.8f, 1.0f));
    if (physics.isOpticallyPumped())
        row("dN_opt", physics.getExcessCarrierDensity(), "%.3e cm^-3",
            ImVec4(0.55f, 0.9f, 1.0f, 1.0f));
    row("dN_drift", physics.getDriftDiffusionExcess(), "%.3e cm^-3",
        ImVec4(1.0f, 0.7f, 0.45f, 1.0f));

    ImGui::SeparatorText("Transport");
    row("mu_n",   physics.getElectronMobility(),    "%8.1f cm^2/Vs");
    row("mu_p",   physics.getHoleMobility(),        "%8.1f cm^2/Vs");
    row("sigma",  physics.getConductivity(),        "%.3e S/cm",
        ImVec4(0.55f, 1.0f, 0.7f, 1.0f));
    row("rho",    physics.getResistivity(),         "%.3e Ohm.cm");

    ImGui::SeparatorText("Magnetic / Hall");
    row("B",      physics.getMagneticField(),       "%8.3f T");
    row("R_H",    physics.getHallCoefficient(),     "%.3e cm^3/C",
        ImVec4(1.0f, 0.7f, 0.45f, 1.0f));

    if (physics.getDopingType() != DopingType::Intrinsic) {
        ImGui::SeparatorText("Ionization");
        row("Ion%",  physics.getIonizationFraction() * 100.0,
            "%6.2f %%");
    }

    // ---- Phase 6: thermal + BJT readouts -------------------------------
    ImGui::SeparatorText("Thermal grid");
    row("T_avg",  dd.meanTemperature(),       "%8.2f K",
        ImVec4(0.55f, 0.85f, 1.00f, 1.0f));
    row("T_peak", dd.maxTemperature(),        "%8.2f K",
        ImVec4(1.00f, 0.55f, 0.30f, 1.0f));
    row("dT_avg", dd.deltaTaverage(),         "%+8.2f K");

    if (dd.deviceMode() == DeviceMode::NpnBjt) {
        ImGui::SeparatorText("BJT (NPN)");
        row("V_BE",  dd.vBE(), "%8.3f V");
        row("V_CE",  dd.vCE(), "%8.3f V");
        row("I_C",   dd.collectorCurrent(), "%.3e a.u.",
            ImVec4(0.60f, 1.00f, 0.80f, 1.0f));
    }

    ImGui::End();
}


// =============================================================================
// Live oscilloscope
// =============================================================================
void drawLiveOscilloscope(const PhysicsEngine& physics, float dt) {
    if (!ImGui::Begin("Live Oscilloscope")) { ImGui::End(); return; }

    static std::vector<float> ts, sigmas, ns, ps;
    static float t_now = 0.0f;
    t_now += dt;
    ts    .push_back(t_now);
    sigmas.push_back(static_cast<float>(physics.getConductivity()));
    ns    .push_back(static_cast<float>(physics.getTotalElectronConc()));
    ps    .push_back(static_cast<float>(physics.getTotalHoleConc()));

    if (ts.size() > kOscilloscopeMaxSamples) {
        const std::size_t excess = ts.size() - kOscilloscopeMaxSamples;
        ts    .erase(ts    .begin(), ts    .begin() + excess);
        sigmas.erase(sigmas.begin(), sigmas.begin() + excess);
        ns    .erase(ns    .begin(), ns    .begin() + excess);
        ps    .erase(ps    .begin(), ps    .begin() + excess);
    }

    if (ImPlot::BeginPlot("Conductivity",
                          ImVec2(-1.0f, ImGui::GetContentRegionAvail().y * 0.55f)))
    {
        ImPlot::SetupAxes("t [s]", "sigma [S/cm]",
                          ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit);
        ImPlot::SetupAxisScale(ImAxis_Y1, ImPlotScale_Log10);
        ImPlot::PlotLine("sigma(t)", ts.data(), sigmas.data(),
                         static_cast<int>(ts.size()));
        ImPlot::EndPlot();
    }

    if (ImPlot::BeginPlot("Carrier concentrations",
                          ImVec2(-1.0f, -1.0f)))
    {
        ImPlot::SetupAxes("t [s]", "carriers [cm^-3]",
                          ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit);
        ImPlot::SetupAxisScale(ImAxis_Y1, ImPlotScale_Log10);
        ImPlot::PlotLine("n(t)", ts.data(), ns.data(),
                         static_cast<int>(ts.size()));
        ImPlot::PlotLine("p(t)", ts.data(), ps.data(),
                         static_cast<int>(ts.size()));
        ImPlot::EndPlot();
    }

    ImGui::End();
}


// =============================================================================
// Crystal view  (RenderTexture embedded as ImGui::Image)
// -----------------------------------------------------------------------------
// Click forwards to DriftDiffusion if interactive mode is on. Image is
// always shown at native size to keep the lattice readable; the scroll
// region inside the ImGui window handles overflow.
// =============================================================================
void drawCrystalViewWindow(const CrystalView&    view,
                           DriftDiffusion&       dd,
                           PhysicsEngine&        physics,
                           UIState&              ui)
{
    if (!ImGui::Begin("Crystal View")) { ImGui::End(); return; }

    const auto size = view.renderTexture().getSize();
    ImGui::Image(view.renderTexture(),
                 sf::Vector2f(static_cast<float>(size.x),
                              static_cast<float>(size.y)));

    // Hover tooltip on the image itself instead of an inline help line.
    if (ui.ddInteractive && ImGui::IsItemHovered()) {
        if (ImGui::BeginTooltip()) {
            ImGui::TextUnformatted(
                "Left-click: deposit a Gaussian generation source here.");
            ImGui::EndTooltip();
        }
    }

    // Click-to-add source.
    if (ui.ddInteractive
        && ImGui::IsItemHovered()
        && ImGui::IsMouseClicked(ImGuiMouseButton_Left))
    {
        const ImVec2 mp     = ImGui::GetMousePos();
        const ImVec2 imgMin = ImGui::GetItemRectMin();
        const ImVec2 imgSz  = ImGui::GetItemRectSize();
        const float u = (mp.x - imgMin.x) / imgSz.x;
        const float v = (mp.y - imgMin.y) / imgSz.y;
        dd.addSource(u, v, ui.sourceIntensity, ui.sourceSigma);
        physics.setDriftDiffusionExcess(dd.globalExcess());
    }

    ImGui::End();
}


// =============================================================================
// Spectrum window  (DOS + absorption coefficient)
// -----------------------------------------------------------------------------
//   * g_c(E) ~ sqrt(E - E_c)     (parabolic CB DOS)
//   * g_v(E) ~ sqrt(E_v - E)     (parabolic VB DOS)
//   * alpha(hv) ~ (hv - E_g)^p
//       p = 1/2 for direct gap (GaAs)
//       p = 2   for indirect gap (Si, Ge); approximated as p = 1 here so
//               the curve fits comfortably on the visible y range.
//
//   Static buffers reused across frames (zero-allocation steady state).
// =============================================================================
void drawSpectrumWindow(const PhysicsEngine& physics) {
    if (!ImGui::Begin("Spectrum")) { ImGui::End(); return; }

    static std::vector<float> energies;
    static std::vector<float> dos_c;
    static std::vector<float> dos_v;
    static std::vector<float> alpha;

    constexpr int kSamples = 200;
    if (energies.capacity() < kSamples + 1) {
        energies.reserve(kSamples + 1);
        dos_c   .reserve(kSamples + 1);
        dos_v   .reserve(kSamples + 1);
        alpha   .reserve(kSamples + 1);
    }
    energies.clear();
    dos_c   .clear();
    dos_v   .clear();
    alpha   .clear();

    const auto& mat = physics.getMaterial();
    const double Eg = physics.getBandgap();
    const double Ec = Eg;
    const double Ev = 0.0;
    const float  p  = mat.isDirectBandgap ? 0.5f : 1.0f;

    for (int i = 0; i <= kSamples; ++i) {
        const double E = -0.4 + 2.4 * i / kSamples;   // -0.4 .. 2.0 eV
        energies.push_back(static_cast<float>(E));

        dos_c.push_back(E > Ec ? static_cast<float>(std::sqrt(E - Ec)) : 0.0f);
        dos_v.push_back(E < Ev ? static_cast<float>(std::sqrt(Ev - E)) : 0.0f);
        alpha.push_back(E > Eg
            ? std::pow(static_cast<float>(E - Eg), p) : 0.0f);
    }

    const float plotH = ImGui::GetContentRegionAvail().y * 0.5f;

    if (ImPlot::BeginPlot("##DOS", ImVec2(-1.0f, plotH))) {
        ImPlot::SetupAxes("E [eV]", "g(E) [a.u.]",
                          ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit);
        ImPlot::PlotLine("g_c (CB)", energies.data(), dos_c.data(), kSamples + 1);
        ImPlot::PlotLine("g_v (VB)", energies.data(), dos_v.data(), kSamples + 1);
        ImPlot::EndPlot();
    }

    if (ImPlot::BeginPlot("##absorption", ImVec2(-1.0f, -1.0f))) {
        ImPlot::SetupAxes("hv [eV]", "alpha(hv) [a.u.]",
                          ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit);
        ImPlot::PlotLine("absorption", energies.data(), alpha.data(),
                         kSamples + 1);
        ImPlot::EndPlot();
    }

    ImGui::End();
}


// =============================================================================
// I-V Curve window  (BJT Gummel + output characteristic)
// -----------------------------------------------------------------------------
// In bulk mode this window stays empty (with a tooltip-only explanation);
// in NPN BJT mode it plots:
//   * a Gummel-style I_C(V_BE) sweep using the engine's emitter-injection
//     model and the current V_CE,
//   * the operating-point marker showing where the user's sliders sit.
//
// Like the spectrum window, the buffers are static and reused frame to
// frame. push_back never reallocates after the first call.
// =============================================================================
void drawIVCurveWindow(const DriftDiffusion& dd) {
    if (!ImGui::Begin("I-V Curve")) { ImGui::End(); return; }

    if (dd.deviceMode() != DeviceMode::NpnBjt) {
        ImGui::TextDisabled(
            "Switch Device mode to NPN BJT to see the I-V characteristics.");
        ImGui::End();
        return;
    }

    static std::vector<float> vbe;
    static std::vector<float> ic;

    constexpr int kSamples = 100;
    if (vbe.capacity() < kSamples + 1) {
        vbe.reserve(kSamples + 1);
        ic .reserve(kSamples + 1);
    }
    vbe.clear();
    ic .clear();

    constexpr float kT_eV = 0.02585f;
    const float V_CE = dd.vCE();

    for (int i = 0; i <= kSamples; ++i) {
        const float v = i * 0.01f;            // 0 .. 1.0 V
        const float arg = std::clamp(v / kT_eV, 0.0f, 25.0f);
        const float n_e = 0.6f * std::exp(arg);
        // Same proxy formula as DriftDiffusion::applyBjtBoundaries.
        const float current = n_e * (V_CE * 0.03f);
        vbe.push_back(v);
        ic .push_back(std::max(current, 1.0e-12f));   // floor for log axis
    }

    if (ImPlot::BeginPlot("##gummel", ImVec2(-1.0f, -1.0f))) {
        ImPlot::SetupAxes("V_BE [V]", "I_C [a.u.]",
                          ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit);
        ImPlot::SetupAxisScale(ImAxis_Y1, ImPlotScale_Log10);
        ImPlot::PlotLine("I_C(V_BE)", vbe.data(), ic.data(), kSamples + 1);

        const float opx[1] = { dd.vBE() };
        const float opy[1] = { std::max(dd.collectorCurrent(), 1.0e-12f) };
        ImPlot::PlotScatter("Operating point", opx, opy, 1);

        ImPlot::EndPlot();
    }
    ImGui::End();
}


// =============================================================================
// Band view (RenderTexture embedded as ImGui::Image)
// =============================================================================
void drawBandViewWindow(const BandView& view) {
    if (!ImGui::Begin("Band Diagram")) { ImGui::End(); return; }
    const auto size = view.renderTexture().getSize();
    ImGui::Image(view.renderTexture(),
                 sf::Vector2f(static_cast<float>(size.x),
                              static_cast<float>(size.y)));
    ImGui::End();
}


// =============================================================================
// Drift-Diffusion inspection
// =============================================================================
void drawCrystalInfoWindow(const DriftDiffusion& dd) {
    if (!ImGui::Begin("Crystal Info")) { ImGui::End(); return; }

    ImGui::SeparatorText("Carrier grid (n)");
    HelpMarker(
        "Solves  dn/dt = D nabla^2 n + G(x,y) - n / tau  on a 2D grid via "
        "explicit FTCS with a per-step CFL clamp. Reflective Neumann at the "
        "walls (BJT contacts override this with Dirichlet boundaries).");
    ImGui::Text("Grid:       %d x %d", dd.width(), dd.height());
    ImGui::Text("n_mean:     %.4f  (a.u.)", dd.meanValue());
    ImGui::Text("n_peak:     %.4f  (a.u.)", dd.maxValue());
    ImGui::Text("Equiv. dN:  %.3e cm^-3",   dd.globalExcess());

    ImGui::SeparatorText("Thermal grid (T)");
    HelpMarker(
        "Solves the heat equation\n"
        "  rho Cp dT/dt = kappa nabla^2 T + H_Joule + H_recomb\n"
        "with Dirichlet edges T = T_ambient (heat-sink contacts). Heat "
        "sources couple electrothermally: more current -> more H_Joule -> "
        "warmer cells -> smaller E_g -> more carriers (the runaway loop).");
    ImGui::Text("T_amb:      %8.2f K", dd.ambientTemperature());
    ImGui::Text("T_mean:     %8.2f K", dd.meanTemperature());
    ImGui::Text("T_peak:     %8.2f K", dd.maxTemperature());
    ImGui::Text("dT_avg:     %+8.2f K", dd.deltaTaverage());

    if (dd.deviceMode() == DeviceMode::NpnBjt) {
        ImGui::SeparatorText("BJT (NPN)");
        ImGui::Text("V_BE:       %5.3f V", dd.vBE());
        ImGui::Text("V_CE:       %5.3f V", dd.vCE());
        ImGui::Text("I_C (proxy):%.3e a.u.", dd.collectorCurrent());
    }

    ImGui::End();
}

} // namespace


// =============================================================================
// Entry point
// =============================================================================
int main() {
    sf::RenderWindow window(
        sf::VideoMode({kWindowWidth, kWindowHeight}),
        kWindowTitle,
        sf::Style::Default);
    window.setFramerateLimit(60);

    if (!ImGui::SFML::Init(window)) {
        std::cerr << "[fatal] ImGui::SFML::Init failed\n";
        return 1;
    }
    ImPlot::CreateContext();
    enableImGuiDocking();
    applyLightStyle();

    sf::Font font;
    (void)loadFont(font);

    auto physics = std::make_unique<PhysicsEngine>(material::Kind::Silicon);
    auto dd      = std::make_unique<DriftDiffusion>(60, 40);
    dd->configureForMaterial(physics->getMaterial());
    dd->setAmbientTemperature(static_cast<float>(physics->getTemperature()));

    auto crystal = std::make_unique<CrystalView>(720);
    auto bands   = std::make_unique<BandView>(720, 540);

    crystal->rebuild(*physics);

    UIState   ui;
    sf::Clock clock;
    bool      running             = true;
    bool      requestResetLayout  = false;

    while (running && window.isOpen()) {
        // ---- Event handling -----------------------------------------------
        while (const std::optional<sf::Event> event = window.pollEvent()) {
            ImGui::SFML::ProcessEvent(window, *event);

            if (event->is<sf::Event::Closed>()) {
                running = false;
            }
            else if (const auto* k = event->getIf<sf::Event::KeyPressed>()) {
                if (k->code == sf::Keyboard::Key::Escape) running = false;
            }
        }

        const sf::Time delta = clock.restart();
        const float    dt    = delta.asSeconds();

        // ---- Sync ambient T into dd (thermal Dirichlet edges) -------------
        dd->setAmbientTemperature(
            static_cast<float>(physics->getTemperature()));

        // ---- Simulation step ----------------------------------------------
        crystal->update(dt, *physics);
        dd->step(dt);
        physics->setDriftDiffusionExcess(dd->globalExcess());

        // Re-roll dopants/carriers periodically when the doping state could
        // have changed (cheap; only resamples carrier visualization).
        // Done implicitly by Controls window callbacks via rebuild path.
        // (Here we just refresh carrier population every frame so optical
        //  and drift-diffusion changes are visible.)
        if (static int frame = 0; (++frame % 30) == 0)
            crystal->rebuild(*physics);

        // ---- Render off-screen targets ------------------------------------
        crystal->render(*physics, *dd, ui.heatmapMode, ui.showVectorField);
        bands->render(*physics, font);

        // ---- ImGui frame --------------------------------------------------
        ImGui::SFML::Update(window, delta);

        (void)beginDockspaceHost(requestResetLayout);
        drawMenuBar(running, ui, *physics, *dd, requestResetLayout);
        ImGui::End();   // close DockHost

        drawControlsWindow   (*physics, *dd, ui);
        drawReadoutsWindow   (*physics, *dd);
        drawLiveOscilloscope (*physics, dt);
        drawSpectrumWindow   (*physics);
        drawIVCurveWindow    (*dd);
        drawCrystalViewWindow(*crystal, *dd, *physics, ui);
        drawBandViewWindow   (*bands);
        drawCrystalInfoWindow(*dd);

        // ---- Final composite ---------------------------------------------
        window.clear(palette::WindowBg);
        ImGui::SFML::Render(window);
        window.display();
    }

    ImPlot::DestroyContext();
    ImGui::SFML::Shutdown();
    return 0;
}
