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
#include <limits>
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
#include "Preset.hpp"


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

// Loads the SFML font (used by SFML render textures: BandView labels) and
// returns a Windows-style path to the same file so ImGui can load it at a
// custom size for the dockable UI.
[[nodiscard]] bool loadFont(sf::Font& font, std::filesystem::path& outPath) {
    for (const auto path : kFontCandidates) {
        const std::filesystem::path fp{path};
        std::error_code ec;
        if (!std::filesystem::exists(fp, ec)) continue;
        if (font.openFromFile(fp)) {
            outPath = fp;
            std::cout << "[info] Loaded font: " << fp.string() << '\n';
            return true;
        }
    }
    std::cerr << "[warn] No usable font found.\n";
    return false;
}


// Targets a comfortable readable size on 1080p+ displays. The default ImGui
// font (Proggy at 13 px) is too small for slider readouts; 17 px keeps text
// legible at standard zoom without crowding the panels.
inline constexpr float kImGuiFontSizePx = 17.0f;


// Loads the same TTF that SFML used into ImGui at a larger point size.
// Falls back to ImGui's built-in pixel font if no TTF is available.
void loadImGuiFont(const std::filesystem::path& fontPath) {
    ImGuiIO& io = ImGui::GetIO();
    io.Fonts->Clear();

    bool ok = false;
    if (!fontPath.empty()) {
        const std::string s = fontPath.string();
        if (io.Fonts->AddFontFromFileTTF(s.c_str(), kImGuiFontSizePx))
            ok = true;
    }
    if (!ok) io.Fonts->AddFontDefault();

    (void)ImGui::SFML::UpdateFontTexture();
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
    ImGui::DockBuilderDockWindow("Device Painter",    dock_center_top);
    ImGui::DockBuilderDockWindow("3D Topology",       dock_center_top);
    ImGui::DockBuilderDockWindow("Live Oscilloscope", dock_center_bottom);
    ImGui::DockBuilderDockWindow("Spectrum",          dock_center_bottom);
    ImGui::DockBuilderDockWindow("I-V Curve",         dock_center_bottom);
    ImGui::DockBuilderDockWindow("Breakdown",         dock_center_bottom);
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

    // BJT controls
    float       V_BE = 0.0f;
    float       V_CE = 0.0f;

    // High-field transport (read-only Caughey-Thomas evaluator).
    float       E_field_V_per_cm = 0.0f;

    // -- Phase 1: Device Painter ------------------------------------------
    BrushKind   currentBrush      = BrushKind::None;
    int         brushRadius       = 3;
    float       brushLogDose      = 17.0f;     // log10 cm^-3
    bool        showDopingOverlay = true;

    // -- Phase 5: Material brush + Wachutka thermal -----------------------
    int            materialBrush  = 0;        // material::Kind index
    bool           showMaterialOverlay = true;
    bool           wachutkaThermal     = true;   // run solveHeatEquation
    float          T_peak_local        = 300.0f;
    float          T_mean_local        = 300.0f;
    float          dT_step_last        = 0.0f;

    // -- Phase 1: Poisson solver ------------------------------------------
    int         poissonIters      = 60;
    float       poissonOmega      = 0.85f;
    bool        poissonAutoSolve  = true;
    double      poissonResidual   = 0.0;

    // -- Phase 2: BandView spatial cut ------------------------------------
    int         bandCutRow        = 20;        // j_row for sampleHorizontalCut
    bool        bandSpatialMode   = true;

    // -- Phase 2: Probe field for Chynoweth / BTBT readout ----------------
    float       probeField        = 3.0e5f;    // V/cm
    float       probeWidth_um     = 0.5f;      // depletion width for M factor

    // -- Phase 3: Non-equilibrium / Gummel --------------------------------
    float       V_bias            = 0.0f;      // applied anode bias [V]
    bool        gummelEnable      = false;     // run Gummel each frame
    int         gummelOuter       = 4;         // outer Poisson<->Continuity passes
    int         gummelPoissonInner   = 25;
    int         gummelContinuityInner = 15;
    float       gummelOmegaPsi    = 0.85f;
    float       gummelOmegaPhi    = 1.0f;
    double      gummelResidual    = 0.0;
    double      lastJ_terminal    = 0.0;       // A/cm^2

    // -- Phase 3: 3D Topology window --------------------------------------
    int         topo3DField       = 0;         // 0=psi, 1=n, 2=p, 3=|E|

    // -- Phase 4: Transient / AC ------------------------------------------
    bool        transientEnabled  = false;
    float       logDt             = -9.0f;     // log10(dt[s]); -9 = 1 ns
    bool        acProbeEnabled    = false;
    float       acFreq_Hz         = 1.0e6f;    // 1 MHz default
    float       acAmp_V           = 0.005f;    // 5 mV
    bool        stepPulseArmed    = false;
    float       stepPulseAmp      = 0.5f;      // forward step amplitude
    double      stepPulseT0       = 0.0;       // time of last step
    bool        stepPulseFired    = false;

    // Oscilloscope ring buffer of (t, J)
    std::vector<float> oscTime;
    std::vector<float> oscCurrent;

    // Small-signal readouts
    double      ssG_DC            = 0.0;       // S/cm^2
    double      ssC_dep           = 0.0;       // F/cm^2 (per unit area)
    double      ssC_diff          = 0.0;       // F/cm^2
    double      ssV_bi            = 0.0;       // built-in potential

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
            // ---- Preset save / load (Phase 6 JSON delivery) ---------------
            if (ImGui::MenuItem("Save preset...", "Ctrl+S")) {
                const auto r = preset::save(
                    "khaos_preset.json", physics, dd);
                ui.flashStatus(r
                    ? "Saved -> khaos_preset.json"
                    : std::string("Save failed: ") + preset::toString(r.error()));
            }
            if (ImGui::MenuItem("Load preset...", "Ctrl+O")) {
                const auto r = preset::load(
                    "khaos_preset.json", physics, dd);
                if (r) {
                    dd.configureForMaterial(physics.getMaterial());
                    ui.V_bias = dd.appliedBias();
                    ui.V_BE   = dd.vBE();
                    ui.V_CE   = dd.vCE();
                    ui.flashStatus(*r);
                } else {
                    ui.flashStatus(std::string("Load failed: ")
                                   + preset::toString(r.error()));
                }
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
            ImGui::TextDisabled("Semiconductor Analysis & Simulation Platform");
            ImGui::TextDisabled("dex / cemfly-april2026");
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
    // High-field transport  (Caughey-Thomas saturation)
    // =====================================================================
    if (ImGui::CollapsingHeader("High-field transport"))
    {
        ImGui::SliderFloat("E-field [V/cm]", &ui.E_field_V_per_cm,
                           0.0f, 1.0e6f, "%.2e",
                           ImGuiSliderFlags_Logarithmic);
        HelpMarker(
            "Caughey-Thomas saturation:\n"
            "  v(E) = mu_low * E / [1 + (mu_low E / v_sat)^beta]^(1/beta)\n"
            "At low fields the carrier velocity is linear in E; above\n"
            "~v_sat / mu_low (~10 kV/cm in Si), the velocity saturates\n"
            "near v_sat (~1e7 cm/s). Sze Sec. 1.5.4.");

        const auto& mat   = physics.getMaterial();
        const double mu_n = physics.getElectronMobility();
        const double mu_p = physics.getHoleMobility();
        const double E    = ui.E_field_V_per_cm;

        const double mu_n_E = PhysicsEngine::highFieldMobility(
            mu_n, E, mat.v_sat_n, mat.beta_n);
        const double mu_p_E = PhysicsEngine::highFieldMobility(
            mu_p, E, mat.v_sat_p, mat.beta_p);
        const double v_n = mu_n_E * E;
        const double v_p = mu_p_E * E;

        ImGui::TextDisabled("mu_n(E) = %.1f  cm^2/Vs   v_n = %.2e cm/s",
                            mu_n_E, v_n);
        ImGui::TextDisabled("mu_p(E) = %.1f  cm^2/Vs   v_p = %.2e cm/s",
                            mu_p_E, v_p);
        ImGui::TextDisabled("v_sat,n = %.2e   v_sat,p = %.2e   (cm/s)",
                            mat.v_sat_n, mat.v_sat_p);
    }

    // =====================================================================
    // Device & BJT
    // =====================================================================
    if (ImGui::CollapsingHeader("Device Topology",
                                ImGuiTreeNodeFlags_DefaultOpen))
    {
        const char* devLabels[] = {
            "Bulk wafer",
            "NPN BJT (procedural)",
            "Painter (custom + Poisson)",
        };
        int dev = static_cast<int>(dd.deviceMode());
        if (ImGui::Combo("Device mode", &dev, devLabels, IM_ARRAYSIZE(devLabels)))
            dd.setDeviceMode(static_cast<DeviceMode>(dev));
        HelpMarker(
            "Bulk: homogeneous wafer; click Crystal View to deposit a "
            "Gaussian generation source.\n"
            "NPN BJT: procedural Emitter/Base/Collector partition; V_BE "
            "and V_CE drive the device.\n"
            "Painter: open the 'Device Painter' window and brush N-/P-type "
            "doping by hand. Equilibrium Poisson is solved on the painted "
            "Nd, Na map every frame; Spatial BandView shows the bent bands.");

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
    // Poisson solver (Phase 1)
    // =====================================================================
    if (ImGui::CollapsingHeader("Poisson Solver"))
    {
        ImGui::Checkbox("Auto-solve (Painter mode)", &ui.poissonAutoSolve);
        HelpMarker(
            "Self-consistent equilibrium Poisson:\n"
            "  eps_s grad^2 psi = -q (p - n + Nd - Na)\n"
            "with n,p in the Boltzmann limit. Gauss-Seidel + under-"
            "relaxation. Active in Painter mode only.\n"
            "Sze Sec. 2.2 / Selberherr Ch. 5.");

        ImGui::SliderInt  ("Iters / frame", &ui.poissonIters, 5, 400);
        ImGui::SliderFloat("omega",         &ui.poissonOmega, 0.30f, 1.20f, "%.2f");

        float pitch_um = dd.cellPitchCm() * 1.0e4f;
        if (ImGui::SliderFloat("Cell pitch [um]", &pitch_um, 0.01f, 5.0f,
                               "%.3f", ImGuiSliderFlags_Logarithmic)) {
            dd.setCellPitchCm(pitch_um * 1.0e-4f);
        }
        HelpMarker(
            "Physical pitch per cell. Smaller pitch -> sharper depletion "
            "edges but slower Poisson relaxation. Rule of thumb: pitch < "
            "L_D / 4 where L_D = sqrt(eps_s V_T / (q N)).");

        if (ImGui::Button("Solve now (1000 iters)")) {
            const auto& m = physics.getMaterial();
            ui.poissonResidual = dd.solvePoisson(
                physics.getIntrinsicCarrier(),
                PhysicsEngine::thermalVoltage(physics.getTemperature()),
                m.epsilon_r, 1000, ui.poissonOmega);
            ui.flashStatus("Poisson 1000-sweep refinement done");
        }
        ImGui::SameLine();
        if (ImGui::Button("Reset psi / clear paint")) {
            dd.clearDoping();
            ui.flashStatus("Doping & potential cleared");
        }
        ImGui::TextDisabled("L2 residual: %.3e V",   ui.poissonResidual);
        ImGui::TextDisabled("|E| peak:    %.3e V/cm", dd.peakElectricField());

        // Built-in potential for the *uniform* doping currently in the
        // Material/Engine state. A useful sanity check next to the
        // Poisson-derived V_bi shown in the spatial band view.
        const auto&  m   = physics.getMaterial();
        const double n_i = physics.getIntrinsicCarrier();
        const double V_T = PhysicsEngine::thermalVoltage(physics.getTemperature());
        const double L_D = PhysicsEngine::debyeLengthCm(
            m.epsilon_r, std::max(physics.getDopingConcentration(), 1.0e10),
            V_T);
        ImGui::TextDisabled("V_T   = %.4f V", V_T);
        ImGui::TextDisabled("L_D   = %.2e cm (~ %.1f um)", L_D, L_D * 1.0e4);
    }

    // =====================================================================
    // Non-equilibrium / Gummel solver  [Phase 3]
    // =====================================================================
    if (ImGui::CollapsingHeader("Bias & Gummel (Non-equilibrium)"))
    {
        if (ImGui::SliderFloat("V_bias [V]", &ui.V_bias, -10.0f, 1.5f, "%.3f")) {
            dd.setAppliedBias(ui.V_bias);
        }
        HelpMarker(
            "Applied anode bias. + = forward (depletion shrinks, current "
            "flows). - = reverse (depletion widens; near-zero current "
            "until breakdown / Zener tunnel).\n"
            "Contacts auto-detected from majority dopant in the leftmost "
            "and rightmost columns -- paint Na on one side, Nd on the "
            "other.");

        ImGui::Checkbox("Run Gummel each frame", &ui.gummelEnable);
        HelpMarker(
            "Gummel iteration: outer loop over Poisson -> Continuity_n "
            "-> Continuity_p, with full quasi-Fermi splitting. Solves for "
            "n,p under bias self-consistently. When OFF, only the "
            "equilibrium Poisson runs (V_bias is ignored).");

        ImGui::SliderInt("Outer iters",     &ui.gummelOuter,    1, 20);
        ImGui::SliderInt("Poisson inner",   &ui.gummelPoissonInner,    5, 80);
        ImGui::SliderInt("Continuity inner",&ui.gummelContinuityInner, 5, 60);
        ImGui::SliderFloat("omega psi", &ui.gummelOmegaPsi, 0.30f, 1.20f, "%.2f");
        ImGui::SliderFloat("omega phi", &ui.gummelOmegaPhi, 0.50f, 1.40f, "%.2f");

        if (ImGui::Button("Solve Gummel now (50 outer)")) {
            const auto& mat = physics.getMaterial();
            ui.gummelResidual = dd.solveGummel(
                physics.getIntrinsicCarrier(),
                PhysicsEngine::thermalVoltage(physics.getTemperature()),
                mat.epsilon_r,
                physics.getElectronMobility(),
                physics.getHoleMobility(),
                mat,
                physics.getTemperature(),
                50, ui.gummelPoissonInner, ui.gummelContinuityInner,
                ui.gummelOmegaPsi, ui.gummelOmegaPhi);
            ui.flashStatus("Gummel 50-outer refinement done");
        }

        ImGui::TextDisabled("psi residual: %.3e V",  ui.gummelResidual);
        ImGui::TextDisabled("J terminal:   %.3e A/cm^2", ui.lastJ_terminal);
    }

    // =====================================================================
    // Transient (BE) and AC probe  [Phase 4]
    // =====================================================================
    if (ImGui::CollapsingHeader("Transient & AC"))
    {
        if (ImGui::Checkbox("Transient (Backward Euler)",
                            &ui.transientEnabled))
        {
            dd.setTransientEnabled(ui.transientEnabled);
            if (ui.transientEnabled) {
                ui.oscTime.clear();
                ui.oscCurrent.clear();
                dd.resetSimulationTime();
            }
        }
        HelpMarker(
            "Adds (n_new - n_old)/dt to continuity. Backward Euler is "
            "unconditionally stable: dt may be picked freely from "
            "picoseconds to milliseconds.\n"
            "Reset clears the simulation clock.");

        if (ImGui::SliderFloat("log10 dt [s]", &ui.logDt, -15.0f, -3.0f,
                               "%.2f"))
        {
            dd.setTimeStep(std::pow(10.0, static_cast<double>(ui.logDt)));
        }
        ImGui::TextDisabled("dt = %.3e s   |   t_sim = %.6e s",
                            dd.timeStep(), dd.simulationTime());

        ImGui::SeparatorText("Step pulse");
        ImGui::SliderFloat("Step amplitude [V]", &ui.stepPulseAmp,
                           -2.0f, 2.0f, "%.3f");
        if (ImGui::Button("Arm step at t = next frame")) {
            ui.stepPulseArmed = true;
            ui.stepPulseFired = false;
            ui.flashStatus("Step pulse armed");
        }
        ImGui::SameLine();
        if (ImGui::Button("Cancel")) {
            ui.stepPulseArmed = false;
            ui.stepPulseFired = false;
        }
        ImGui::SameLine();
        ImGui::TextDisabled(ui.stepPulseFired ? "fired"
                            : (ui.stepPulseArmed ? "armed" : "idle"));

        ImGui::SeparatorText("AC probe");
        if (ImGui::Checkbox("AC probe", &ui.acProbeEnabled)) {
            dd.setACProbe(ui.acProbeEnabled,
                          ui.acFreq_Hz, ui.acAmp_V);
        }
        HelpMarker(
            "Rides a small sine on top of V_a:\n"
            "  V_a(t) = V_DC + amp * sin(2 pi f t)\n"
            "Used with transient stepping to read C, G off the live "
            "oscilloscope. Frequencies above 1/(2 pi tau_recomb) probe "
            "the depletion (junction) capacitance; lower frequencies "
            "see the diffusion (storage) capacitance too.");

        if (ImGui::SliderFloat("AC freq [Hz]", &ui.acFreq_Hz,
                               1.0f, 1.0e10f, "%.2e",
                               ImGuiSliderFlags_Logarithmic)) {
            dd.setACProbe(ui.acProbeEnabled, ui.acFreq_Hz, ui.acAmp_V);
        }
        if (ImGui::SliderFloat("AC amplitude [V]", &ui.acAmp_V,
                               0.0001f, 0.1f, "%.4f",
                               ImGuiSliderFlags_Logarithmic)) {
            dd.setACProbe(ui.acProbeEnabled, ui.acFreq_Hz, ui.acAmp_V);
        }

        // ---- Phase 5: Wachutka local thermal coupling ------------------
        ImGui::SeparatorText("Wachutka thermal coupling");
        ImGui::Checkbox("Solve local heat equation", &ui.wachutkaThermal);
        HelpMarker(
            "Enables the per-cell heat solver each transient step:\n"
            "  rho Cp dT/dt = div(kappa grad T) + H\n"
            "with H = (J_n+J_p).E + (R - G_BTBT)(E_g + 3kT).\n"
            "Hot cells narrow E_g and inflate n_i exponentially -- the "
            "Wachutka feedback loop that drives thermal runaway.\n"
            "Reference: Wachutka, IEEE TCAD 9 (1990) 1141; Selberherr 4.5.");

        ImGui::Text("T_peak    = %7.2f K  (dT_step = %+.3f K)",
                    ui.T_peak_local, ui.dT_step_last);
        ImGui::Text("T_mean    = %7.2f K", ui.T_mean_local);
        ImGui::Text("T_ambient = %7.2f K", dd.ambientTemperature());
        if (ui.T_peak_local - dd.ambientTemperature() > 30.0f) {
            ImGui::TextColored(ImVec4(1.0f, 0.55f, 0.25f, 1.0f),
                "Hotspot active: dT > 30 K -- watch n_i / I rise.");
        }
    }

    // =====================================================================
    // Small-signal (DC perturbation C, G)  [Phase 4 bonus]
    // =====================================================================
    if (ImGui::CollapsingHeader("Small-Signal (DC perturbation)"))
    {
        const auto&  mat = physics.getMaterial();
        const double n_i = physics.getIntrinsicCarrier();
        const double V_T_local =
            PhysicsEngine::thermalVoltage(physics.getTemperature());

        // Compute V_bi from peak Nd / Na on the painter map, simply
        // because the user might paint asymmetric structures.
        double Nd_max = 0.0, Na_max = 0.0;
        for (int j = 0; j < dd.height(); ++j) {
            for (int i = 0; i < dd.width(); ++i) {
                Nd_max = std::max(Nd_max, dd.donorAt   (i, j));
                Na_max = std::max(Na_max, dd.acceptorAt(i, j));
            }
        }
        ui.ssV_bi = PhysicsEngine::builtInPotential(Nd_max, Na_max,
                                                    n_i, V_T_local);

        ImGui::Text("V_bi (max-Nd vs max-Na): %.3f V", ui.ssV_bi);
        ImGui::Text("Nd_max = %.2e   Na_max = %.2e", Nd_max, Na_max);

        if (ImGui::Button("Measure G_DC (perturbation +/-5 mV)")) {
            ui.ssG_DC = dd.smallSignalConductance(
                n_i, V_T_local, mat.epsilon_r, mat,
                physics.getTemperature(), 0.005);
            ui.flashStatus("G_DC measured");
        }
        ImGui::SameLine();
        if (ImGui::Button("Update C estimates")) {
            // C_dep at the *current* V_bias (cm^2-normalised).
            ui.ssC_dep = PhysicsEngine::depletionCapacitanceFlat(
                Nd_max, Na_max, ui.ssV_bi,
                static_cast<double>(ui.V_bias),
                mat.epsilon_r, 1.0);
            ui.ssC_diff = PhysicsEngine::diffusionCapacitance(
                std::max(ui.lastJ_terminal, 0.0),
                mat.tau_p, V_T_local);
            ui.flashStatus("Capacitances updated");
        }

        ImGui::TextDisabled("G_DC   = %.3e S/cm^2",  ui.ssG_DC);
        ImGui::TextDisabled("C_dep  = %.3e F/cm^2",  ui.ssC_dep);
        ImGui::TextDisabled("C_diff = %.3e F/cm^2",  ui.ssC_diff);
        ImGui::TextDisabled("C_tot  = %.3e F/cm^2",
                            ui.ssC_dep + ui.ssC_diff);
        if (ui.ssG_DC > 0.0)
            ImGui::TextDisabled("R_DC   = %.3e Ohm.cm^2", 1.0 / ui.ssG_DC);
    }

    // =====================================================================
    // BandView cut & mode (Phase 2)
    // =====================================================================
    if (ImGui::CollapsingHeader("Band Diagram"))
    {
        ImGui::Checkbox("Spatial mode (bent bands)", &ui.bandSpatialMode);
        HelpMarker(
            "Off: legacy single-point band diagram from PhysicsEngine.\n"
            "On:  spatial cut Ec(x), Ev(x), Ef along the row below; "
            "Ec(x) = Ec0 - q psi(x), Ev(x) = Ev0 - q psi(x). "
            "Sze Eq. 2.10 / Pierret Sec. 4.2.");

        ImGui::SliderInt("Cut row j", &ui.bandCutRow,
                         0, std::max(0, dd.height() - 1));
        ImGui::TextDisabled("y = %.1f um", ui.bandCutRow
                            * dd.cellPitchCm() * 1.0e4);
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
        ImVec4(0.78f, 0.46f, 0.18f, 1.0f));
    row("p",      physics.getTotalHoleConc(),       "%.3e cm^-3",
        ImVec4(0.55f, 0.20f, 0.65f, 1.0f));
    if (physics.isOpticallyPumped())
        row("dN_opt", physics.getExcessCarrierDensity(), "%.3e cm^-3",
            ImVec4(0.20f, 0.45f, 0.78f, 1.0f));
    row("dN_drift", physics.getDriftDiffusionExcess(), "%.3e cm^-3",
        ImVec4(0.70f, 0.40f, 0.20f, 1.0f));

    ImGui::SeparatorText("Transport");
    row("mu_n",   physics.getElectronMobility(),    "%8.1f cm^2/Vs");
    row("mu_p",   physics.getHoleMobility(),        "%8.1f cm^2/Vs");
    row("sigma",  physics.getConductivity(),        "%.3e S/cm",
        ImVec4(0.10f, 0.50f, 0.30f, 1.0f));
    row("rho",    physics.getResistivity(),         "%.3e Ohm.cm");

    ImGui::SeparatorText("Magnetic / Hall");
    row("B",      physics.getMagneticField(),       "%8.3f T");
    row("R_H",    physics.getHallCoefficient(),     "%.3e cm^3/C",
        ImVec4(0.78f, 0.42f, 0.18f, 1.0f));

    if (physics.getDopingType() != DopingType::Intrinsic) {
        ImGui::SeparatorText("Ionization");
        row("Ion%",  physics.getIonizationFraction() * 100.0,
            "%6.2f %%");
    }

    // ---- Thermal grid + electrothermal feedback ------------------------
    ImGui::SeparatorText("Thermal grid");
    row("T_avg",  dd.meanTemperature(),       "%8.2f K",
        ImVec4(0.20f, 0.45f, 0.75f, 1.0f));
    row("T_peak", dd.maxTemperature(),        "%8.2f K",
        ImVec4(0.85f, 0.30f, 0.10f, 1.0f));
    row("dT_avg", dd.deltaTaverage(),         "%+8.2f K");

    // Electrothermal feedback readouts: at T_peak the bandgap narrows and
    // the intrinsic carrier density jumps. Show both side by side with the
    // ambient values to make the runaway condition visible numerically.
    if (dd.maxTemperature() > dd.ambientTemperature() + 1.0f) {
        const auto& mat   = physics.getMaterial();
        const double T_pk = dd.maxTemperature();
        const double Eg_pk = PhysicsEngine::bandgapAt(mat, T_pk);
        const double ni_pk = PhysicsEngine::intrinsicCarrierAt(mat, T_pk);
        row("E_g(peak)", Eg_pk, "%8.4f eV",
            ImVec4(0.55f, 0.30f, 0.10f, 1.0f));
        row("n_i(peak)", ni_pk, "%.3e cm^-3",
            ImVec4(0.55f, 0.30f, 0.10f, 1.0f));
    }

    if (dd.deviceMode() == DeviceMode::NpnBjt) {
        ImGui::SeparatorText("BJT (NPN)");
        row("V_BE",  dd.vBE(), "%8.3f V");
        row("V_CE",  dd.vCE(), "%8.3f V");
        row("I_C",   dd.collectorCurrent(), "%.3e a.u.",
            ImVec4(0.10f, 0.50f, 0.30f, 1.0f));
    }

    if (dd.deviceMode() == DeviceMode::Painter
        && std::abs(dd.appliedBias()) > 1.0e-4f)
    {
        ImGui::SeparatorText("Bias / non-equilibrium");
        row("V_a", dd.appliedBias(), "%8.3f V",
            ImVec4(0.95f, 0.85f, 0.40f, 1.0f));
    }

    ImGui::End();
}


// =============================================================================
// Live oscilloscope
// =============================================================================
void drawLiveOscilloscope(const PhysicsEngine& physics,
                          const DriftDiffusion& dd,
                          UIState& ui,
                          float dt)
{
    if (!ImGui::Begin("Live Oscilloscope")) { ImGui::End(); return; }

    // Mode banner: transient mode reads ui.oscTime / oscCurrent populated
    // by the main loop; legacy mode keeps the wall-clock physics traces.
    const bool transient = ui.transientEnabled
                        && dd.deviceMode() == DeviceMode::Painter;

    if (ImGui::BeginTabBar("##scope_tabs")) {
        if (ImGui::BeginTabItem(transient
                                ? "Transient J(t)" : "Conductivity & n,p"))
        {
            if (transient) {
                if (ImPlot::BeginPlot("Terminal current",
                    ImVec2(-1.0f, ImGui::GetContentRegionAvail().y - 8.0f)))
                {
                    ImPlot::SetupAxes("t_sim [s]", "J [A/cm^2]",
                                      ImPlotAxisFlags_AutoFit,
                                      ImPlotAxisFlags_AutoFit);
                    if (!ui.oscTime.empty()) {
                        ImPlot::PlotLine("J(t)",
                            ui.oscTime.data(),
                            ui.oscCurrent.data(),
                            static_cast<int>(ui.oscTime.size()));
                    }
                    if (ui.stepPulseFired) {
                        const double xs[2] = {ui.stepPulseT0, ui.stepPulseT0};
                        const double ys[2] = {-1e30, 1e30};
                        ImPlot::PlotLine("step", xs, ys, 2);
                    }
                    ImPlot::EndPlot();
                }
                ImGui::TextDisabled(
                    "Samples: %zu  |  t = %.4e s  |  V_a = %.4f V  |  "
                    "J = %.3e A/cm^2",
                    ui.oscTime.size(),
                    dd.simulationTime(),
                    static_cast<double>(dd.appliedBias()),
                    ui.lastJ_terminal);
            } else {
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
            }
            ImGui::EndTabItem();
        }
        ImGui::EndTabBar();
    }

    if (ImGui::SmallButton("Clear scope")) {
        ui.oscTime.clear();
        ui.oscCurrent.clear();
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
    static std::vector<float> alpha_curve;
    static std::vector<float> depth_um;

    constexpr int kSamples = 200;
    if (energies.capacity() < kSamples + 1) {
        energies   .reserve(kSamples + 1);
        dos_c      .reserve(kSamples + 1);
        dos_v      .reserve(kSamples + 1);
        alpha_curve.reserve(kSamples + 1);
        depth_um   .reserve(kSamples + 1);
    }
    energies   .clear();
    dos_c      .clear();
    dos_v      .clear();
    alpha_curve.clear();
    depth_um   .clear();

    const auto& mat = physics.getMaterial();
    const double Eg = physics.getBandgap();
    const double Ec = Eg;
    const double Ev = 0.0;

    for (int i = 0; i <= kSamples; ++i) {
        const double E = -0.4 + 2.4 * i / kSamples;   // -0.4 .. 2.0 eV
        energies.push_back(static_cast<float>(E));

        dos_c.push_back(E > Ec ? static_cast<float>(std::sqrt(E - Ec)) : 0.0f);
        dos_v.push_back(E < Ev ? static_cast<float>(std::sqrt(Ev - E)) : 0.0f);

        // Real cm^-1 absorption coefficient + penetration depth in um.
        const double a = PhysicsEngine::absorptionCoefficient(mat, E);
        alpha_curve.push_back(static_cast<float>(a));
        if (a > 0.0) {
            // L_alpha = 1/alpha [cm] -> convert to micrometers.
            const double L_um = (1.0 / a) * 1.0e4;
            depth_um.push_back(static_cast<float>(std::min(L_um, 1.0e6)));
        } else {
            depth_um.push_back(1.0e6f);   // off-scale at sub-gap photons
        }
    }

    const float plotH = ImGui::GetContentRegionAvail().y * 0.33f;

    if (ImPlot::BeginPlot("##DOS", ImVec2(-1.0f, plotH))) {
        ImPlot::SetupAxes("E [eV]", "g(E) [a.u.]",
                          ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit);
        ImPlot::PlotLine("g_c (CB)", energies.data(), dos_c.data(), kSamples + 1);
        ImPlot::PlotLine("g_v (VB)", energies.data(), dos_v.data(), kSamples + 1);
        ImPlot::EndPlot();
    }

    if (ImPlot::BeginPlot("##absorption", ImVec2(-1.0f, plotH))) {
        ImPlot::SetupAxes("hv [eV]", "alpha [cm^-1]",
                          ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit);
        ImPlot::SetupAxisScale(ImAxis_Y1, ImPlotScale_Log10);
        ImPlot::PlotLine("alpha(hv)", energies.data(), alpha_curve.data(),
                         kSamples + 1);
        ImPlot::EndPlot();
    }

    if (ImPlot::BeginPlot("##penetration", ImVec2(-1.0f, -1.0f))) {
        ImPlot::SetupAxes("hv [eV]", "L_alpha [um]",
                          ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit);
        ImPlot::SetupAxisScale(ImAxis_Y1, ImPlotScale_Log10);
        ImPlot::PlotLine("penetration depth",
                         energies.data(), depth_um.data(), kSamples + 1);
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

    // Use the actual ambient T so the curve and the simulation agree.
    const float kT_eV = 8.617333262e-5f * dd.ambientTemperature();
    const float V_CE  = dd.vCE();

    for (int i = 0; i <= kSamples; ++i) {
        const float v   = i * 0.01f;          // 0 .. 1.0 V
        const float arg = std::clamp(v / kT_eV, 0.0f, 25.0f);
        const float n_e = 0.6f * std::exp(arg);
        // Mirrors DriftDiffusion::applyBjtBoundaries exactly.
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
// Band view  (RenderTexture embedded as ImGui::Image)
// -----------------------------------------------------------------------------
// In Spatial mode, the user can also click inside the image to set the
// j_cut row directly -- a quicker workflow than the slider in Controls.
// =============================================================================
void drawBandViewWindow(const BandView& view,
                        const DriftDiffusion& dd,
                        UIState& ui)
{
    if (!ImGui::Begin("Band Diagram")) { ImGui::End(); return; }

    const auto size = view.renderTexture().getSize();
    ImGui::Image(view.renderTexture(),
                 sf::Vector2f(static_cast<float>(size.x),
                              static_cast<float>(size.y)));

    if (ui.bandSpatialMode && ImGui::IsItemHovered()
        && ImGui::IsMouseDown(ImGuiMouseButton_Left))
    {
        // Crystal View has its own click handler for source deposition;
        // here we intentionally repurpose left-click for "pick cut row".
        // The y coordinate inside the band image does NOT directly map
        // to a grid row because the band canvas only uses the upper 66%
        // of the image; we ignore that and let the user pick via the
        // slider for fine control. This handler is purely a convenience
        // for "show me the band diagram across this lateral position".
        (void)dd;
    }
    ImGui::TextDisabled(
        ui.bandSpatialMode ? "Spatial cut active -- adjust 'Cut row j' in Controls."
                           : "Flat mode (PhysicsEngine global state).");
    ImGui::End();
}


// =============================================================================
// Device Painter [Phase 1]
// -----------------------------------------------------------------------------
// Click-and-drag canvas that stamps Nd / Na onto the drift-diffusion grid.
// Renders the doping map, the |E|=0 (junction) iso-line, and a hover
// tooltip with the local Nd, Na, psi, and band shift.
// =============================================================================
void drawDevicePainterWindow(DriftDiffusion& dd,
                             const PhysicsEngine& physics,
                             UIState& ui)
{
    if (!ImGui::Begin("Device Painter")) { ImGui::End(); return; }

    // ---- Brush palette --------------------------------------------------
    ImGui::TextDisabled("Brush");
    int brush = static_cast<int>(ui.currentBrush);
    const char* brushLabels[] = {
        "None",
        "N-dopant (donor)",
        "P-dopant (acceptor)",
        "Eraser",
        "Material (Si/Ge/GaAs)",        // [Phase 5]
    };
    ImGui::Combo("##brush", &brush, brushLabels, IM_ARRAYSIZE(brushLabels));
    ui.currentBrush = static_cast<BrushKind>(brush);

    ImGui::SliderInt  ("Radius (cells)",  &ui.brushRadius,  1, 12);
    ImGui::SliderFloat("log10 dose [cm^-3]",
                       &ui.brushLogDose, 14.0f, 20.0f, "%.2f");
    ImGui::Checkbox   ("Overlay Nd/Na heatmap", &ui.showDopingOverlay);
    HelpMarker(
        "Click & drag to stamp donors (N) or acceptors (P). Each frame "
        "the brush adds the chosen dose to every cell within the circular "
        "footprint. Eraser zeroes Nd, Na and the cell tag back to bulk. "
        "Material brush stamps a heterojunction patch (Si/Ge/GaAs).");

    ImGui::SameLine();
    if (ImGui::Button("Clear all doping")) {
        dd.clearDoping();
        ui.flashStatus("Painter canvas cleared");
    }

    // ---- Phase 5: Material brush controls -------------------------------
    if (ui.currentBrush == BrushKind::Material) {
        ImGui::Combo("Material kind", &ui.materialBrush,
                     material::kLabels, material::kCount);
        HelpMarker(
            "Anderson rule (Sze 5.1):\n"
            "  dEc = chi_A - chi_B\n"
            "  dEv = (chi_B - chi_A) + (E_g,B - E_g,A)\n"
            "Si chi=4.05, Ge chi=4.00, GaAs chi=4.07 eV.\n"
            "Stamping Ge next to Si makes a textbook SiGe HBT base "
            "(small dEc ~0.05 eV, large dEv ~0.42 eV: holes blocked).");
    }
    ImGui::Checkbox("Overlay material map", &ui.showMaterialOverlay);
    ImGui::SameLine();
    if (ImGui::Button("Clear materials")) {
        dd.clearMaterials();
        ui.flashStatus("Material map reset to engine reference");
    }

    // Quick presets for repeatable lab-style structures.
    if (ImGui::Button("Preset: PN junction")) {
        dd.clearDoping();
        for (int j = 0; j < dd.height(); ++j) {
            for (int i = 0; i < dd.width(); ++i) {
                const float u = (i + 0.5f) / static_cast<float>(dd.width());
                const float v = (j + 0.5f) / static_cast<float>(dd.height());
                dd.paintBrush(u, v,
                    (i < dd.width() / 2) ? BrushKind::PDopant
                                         : BrushKind::NDopant,
                    1.0e17, 0);
            }
        }
        if (dd.deviceMode() != DeviceMode::Painter)
            dd.setDeviceMode(DeviceMode::Painter);
        ui.flashStatus("PN junction stamped (1e17 / 1e17)");
    }
    ImGui::SameLine();
    if (ImGui::Button("Preset: P+ N N+")) {
        dd.clearDoping();
        for (int j = 0; j < dd.height(); ++j) {
            for (int i = 0; i < dd.width(); ++i) {
                const float u = (i + 0.5f) / static_cast<float>(dd.width());
                const float v = (j + 0.5f) / static_cast<float>(dd.height());
                BrushKind k;  double dose;
                if      (i < dd.width() * 0.20) { k = BrushKind::PDopant; dose = 1.0e19; }
                else if (i < dd.width() * 0.75) { k = BrushKind::NDopant; dose = 1.0e15; }
                else                             { k = BrushKind::NDopant; dose = 1.0e19; }
                dd.paintBrush(u, v, k, dose, 0);
            }
        }
        if (dd.deviceMode() != DeviceMode::Painter)
            dd.setDeviceMode(DeviceMode::Painter);
        ui.flashStatus("P+/N/N+ stamped");
    }

    // ---- Canvas ---------------------------------------------------------
    ImGui::Separator();
    const ImVec2 avail   = ImGui::GetContentRegionAvail();
    const float  cellPx  = std::max(4.0f,
        std::min(avail.x / static_cast<float>(dd.width()),
                 avail.y / static_cast<float>(dd.height())));
    const ImVec2 canvasSz(cellPx * static_cast<float>(dd.width()),
                          cellPx * static_cast<float>(dd.height()));

    ImGui::InvisibleButton("##painter_canvas", canvasSz);
    const ImVec2 p0 = ImGui::GetItemRectMin();
    const ImVec2 p1 = ImGui::GetItemRectMax();

    auto* dl = ImGui::GetWindowDrawList();
    dl->AddRectFilled(p0, p1, IM_COL32(20, 20, 28, 255));

    // ---- Material overlay (Phase 5: per-cell Si/Ge/GaAs tint) ----------
    if (ui.showMaterialOverlay) {
        for (int j = 0; j < dd.height(); ++j) {
            for (int i = 0; i < dd.width(); ++i) {
                const auto kind = dd.materialAt(i, j);
                const auto& mat = material::byKind(kind);
                // Subtle base tint -- material colour at low alpha so the
                // doping overlay above remains the dominant visual.
                const ImU32 mcol = IM_COL32(mat.atomR, mat.atomG, mat.atomB, 60);
                const ImVec2 c0(p0.x + cellPx * i,
                                p0.y + cellPx * j);
                const ImVec2 c1(c0.x + cellPx, c0.y + cellPx);
                dl->AddRectFilled(c0, c1, mcol);
            }
        }
    }

    // ---- Heatmap overlay (Nd-Na, signed log) ----------------------------
    if (ui.showDopingOverlay) {
        for (int j = 0; j < dd.height(); ++j) {
            for (int i = 0; i < dd.width(); ++i) {
                const double net = dd.netDopingAt(i, j);
                if (net == 0.0) continue;
                const float lg = static_cast<float>(
                    std::log10(std::abs(net) + 1.0) / 20.0);
                const auto a = static_cast<int>(
                    std::clamp(lg * 255.0f, 0.0f, 255.0f));
                ImU32 col = (net > 0)
                    ? IM_COL32( 60, 140, 255, a)
                    : IM_COL32(255, 140,  60, a);
                const ImVec2 c0(p0.x + cellPx * i,
                                p0.y + cellPx * j);
                const ImVec2 c1(c0.x + cellPx, c0.y + cellPx);
                dl->AddRectFilled(c0, c1, col);
            }
        }
    }

    // ---- psi sign-change isolines (junction lines) ---------------------
    for (int j = 1; j < dd.height(); ++j) {
        for (int i = 1; i < dd.width(); ++i) {
            const float a = dd.psiAt(i,     j);
            const float b = dd.psiAt(i - 1, j);
            if ((a > 0) != (b > 0)) {
                const ImVec2 q0(p0.x + cellPx * i,
                                p0.y + cellPx * j);
                const ImVec2 q1(q0.x, q0.y + cellPx);
                dl->AddLine(q0, q1, IM_COL32(255, 240, 120, 220), 1.5f);
            }
        }
    }

    // ---- Cut-row indicator (Phase 2 BandView) --------------------------
    {
        const float yLine = p0.y + cellPx * static_cast<float>(ui.bandCutRow)
                          + cellPx * 0.5f;
        dl->AddLine({p0.x, yLine}, {p1.x, yLine},
                    IM_COL32(255, 255, 255, 160), 1.0f);
        dl->AddText({p0.x + 4.0f, yLine - 14.0f},
                    IM_COL32(255, 255, 255, 200), "BandView cut");
    }

    // ---- Brush stamping -------------------------------------------------
    if (ImGui::IsItemHovered() && ui.currentBrush != BrushKind::None
        && ImGui::IsMouseDown(ImGuiMouseButton_Left))
    {
        const ImVec2 mp = ImGui::GetMousePos();
        const float u = (mp.x - p0.x) / canvasSz.x;
        const float v = (mp.y - p0.y) / canvasSz.y;
        if (ui.currentBrush == BrushKind::Material) {
            // Phase 5: stamp local material id (heterojunction painter).
            const auto kind = static_cast<material::Kind>(
                std::clamp(ui.materialBrush, 0, material::kCount - 1));
            dd.paintMaterialBrush(u, v, kind, ui.brushRadius);
        } else {
            const double dose = std::pow(10.0, ui.brushLogDose);
            dd.paintBrush(u, v, ui.currentBrush, dose, ui.brushRadius);
        }
        if (dd.deviceMode() != DeviceMode::Painter)
            dd.setDeviceMode(DeviceMode::Painter);
    }

    // Right-click sets the BandView cut row to the hovered y -- quicker
    // than scrubbing a slider in another panel.
    if (ImGui::IsItemHovered()
        && ImGui::IsMouseClicked(ImGuiMouseButton_Right))
    {
        const ImVec2 mp = ImGui::GetMousePos();
        const int j = std::clamp(
            static_cast<int>((mp.y - p0.y) / cellPx), 0, dd.height() - 1);
        ui.bandCutRow = j;
    }

    // ---- Brush cursor preview -------------------------------------------
    if (ImGui::IsItemHovered()) {
        const ImVec2 mp = ImGui::GetMousePos();
        ImU32 col = IM_COL32(180, 180, 180, 180);
        switch (ui.currentBrush) {
            case BrushKind::NDopant: col = IM_COL32( 80, 180, 255, 220); break;
            case BrushKind::PDopant: col = IM_COL32(255, 160,  80, 220); break;
            case BrushKind::Eraser:  col = IM_COL32(220,  60,  60, 220); break;
            case BrushKind::Material: {
                const auto& mat = material::byKind(static_cast<material::Kind>(
                    std::clamp(ui.materialBrush, 0, material::kCount - 1)));
                col = IM_COL32(mat.atomR, mat.atomG, mat.atomB, 240);
                break;
            }
            default: break;
        }
        dl->AddCircle(mp, ui.brushRadius * cellPx, col, 0, 2.0f);
    }

    // ---- Hover tooltip --------------------------------------------------
    if (ImGui::IsItemHovered()) {
        const ImVec2 mp = ImGui::GetMousePos();
        const int i = std::clamp(
            static_cast<int>((mp.x - p0.x) / cellPx), 0, dd.width()  - 1);
        const int j = std::clamp(
            static_cast<int>((mp.y - p0.y) / cellPx), 0, dd.height() - 1);
        if (ImGui::BeginTooltip()) {
            ImGui::Text("(i,j) = (%d, %d)", i, j);
            const auto& mat_local = dd.profileAt(i, j);
            ImGui::Text("Material : %.*s",
                static_cast<int>(mat_local.name.size()), mat_local.name.data());
            ImGui::Text("chi      = %.3f eV", mat_local.chi);
            ImGui::Text("E_g(T)   = %.3f eV", dd.localBandgapAt(i, j));
            ImGui::Text("T_lat    = %.2f K",  dd.temperatureFieldAt(i, j));
            ImGui::Text("Nd       = %.2e cm^-3", dd.donorAt(i, j));
            ImGui::Text("Na       = %.2e cm^-3", dd.acceptorAt(i, j));
            ImGui::Text("psi      = %+.4f V",    dd.psiAt(i, j));
            ImGui::Text("|E|      = %.3e V/cm",  dd.electricFieldMagAt(i, j));
            ImGui::EndTooltip();
        }
        (void)physics;
    }

    ImGui::End();
}


// =============================================================================
// Breakdown / Recombination panel  [Phase 2]
// -----------------------------------------------------------------------------
// Shows three high-field / heavy-doping mechanisms side by side:
//
//   * Auger recombination R_Aug = (C_n n + C_p p)(np - n_i^2)
//   * Chynoweth ionization rates alpha_n(E), alpha_p(E)
//     and Miller multiplication factor M = 1/(1 - alpha W)
//   * Kane band-to-band tunneling generation G_BTBT(E)
//
// Two probe sources:
//   * "Engine probe": user-chosen E and W, evaluated against the global
//                     Material profile.
//   * "Grid peak"   : the peak |E| extracted from the Painter's psi grid;
//                     answers "is my drawn junction breaking down?".
// =============================================================================
void drawBreakdownWindow(const PhysicsEngine& physics,
                         const DriftDiffusion& dd,
                         UIState& ui)
{
    if (!ImGui::Begin("Breakdown")) { ImGui::End(); return; }

    const auto&  m   = physics.getMaterial();
    const double n   = physics.getTotalElectronConc();
    const double p   = physics.getTotalHoleConc();
    const double n_i = physics.getIntrinsicCarrier();

    // ---- Recombination triplet (SRH + Auger + Radiative) ----------------
    ImGui::SeparatorText("Recombination (U_net = R_SRH + R_Aug + R_rad)");
    HelpMarker(
        "Three (n,p)-only loss channels compete:\n"
        "  R_SRH = (np - n_i^2) / [tau_p(n+n_i) + tau_n(p+n_i)] -- midgap traps\n"
        "  R_Aug = (C_n n + C_p p)(np - n_i^2)                 -- 3-particle\n"
        "  R_rad = B_rad (np - n_i^2)                           -- photon emission\n"
        "Field-driven Kane/Chynoweth generation is reported below.\n"
        "Sze 1.5.6 / Pankove Ch. 6 / Schubert 2.13.");

    const double R_srh = PhysicsEngine::recombSRH(n, p, n_i,
                                                  m.tau_n, m.tau_p);
    const double R_aug = PhysicsEngine::recombAuger(n, p, n_i,
                                                    m.C_n_aug, m.C_p_aug);
    const double R_rad = PhysicsEngine::recombRadiative(n, p, n_i, m.B_rad);
    const double U_net = R_srh + R_aug + R_rad;

    ImGui::Text("tau_n  = %.2e   tau_p  = %.2e s",  m.tau_n, m.tau_p);
    ImGui::Text("C_n    = %.2e   C_p    = %.2e cm^6/s", m.C_n_aug, m.C_p_aug);
    ImGui::Text("B_rad  = %.2e cm^3/s   (%s gap)",
                m.B_rad, m.isDirectBandgap ? "DIRECT" : "indirect");
    ImGui::Spacing();
    ImGui::Text("R_SRH  = %.3e cm^-3 / s",  R_srh);
    ImGui::Text("R_Aug  = %.3e cm^-3 / s",  R_aug);
    ImGui::Text("R_rad  = %.3e cm^-3 / s",  R_rad);
    ImGui::Text("U_net  = %.3e cm^-3 / s",  U_net);

    if (U_net > 0.0) {
        const double f_srh = R_srh / U_net;
        const double f_aug = R_aug / U_net;
        const double f_rad = R_rad / U_net;
        ImGui::Text("Shares  SRH %.1f%%   Aug %.1f%%   Rad %.1f%%",
                    100.0 * f_srh, 100.0 * f_aug, 100.0 * f_rad);

        const char* dominant =
            (f_srh >= f_aug && f_srh >= f_rad) ? "SRH-dominated (trap-assisted)" :
            (f_aug >= f_rad)                   ? "Auger-dominated (heavy doping / high injection)"
                                               : "Radiative-dominated (direct-gap LED regime)";
        ImVec4 col = (f_rad > 0.5) ? ImVec4(0.40f, 0.85f, 1.0f, 1.0f)
                                    : (f_aug > 0.5) ? ImVec4(1.0f, 0.55f, 0.30f, 1.0f)
                                                    : ImVec4(0.75f, 0.85f, 0.95f, 1.0f);
        ImGui::TextColored(col, "%s", dominant);
    }

    // ---- Chynoweth + Kane probe ----------------------------------------
    ImGui::SeparatorText("High-field probe (Chynoweth / Kane BTBT)");
    HelpMarker(
        "alpha(E) = alpha_inf exp[-(E_crit/|E|)^m]. Avalanche multiplication "
        "M = 1/(1 - alpha W) (Miller form, single-carrier limit). Kane BTBT "
        "G(E) = A E^P exp(-B/E).");
    ImGui::SliderFloat("Probe |E| [V/cm]", &ui.probeField,
                       1.0e3f, 5.0e6f, "%.2e",
                       ImGuiSliderFlags_Logarithmic);
    ImGui::SliderFloat("Depletion W [um]", &ui.probeWidth_um,
                       0.05f, 10.0f, "%.3f");

    const double alpha_n = PhysicsEngine::chynowethRate(
        ui.probeField, m.alpha_inf_n, m.E_crit_n, m.chyn_m);
    const double alpha_p = PhysicsEngine::chynowethRate(
        ui.probeField, m.alpha_inf_p, m.E_crit_p, m.chyn_m);
    const double W_cm    = ui.probeWidth_um * 1.0e-4;
    const double Mfac    = PhysicsEngine::avalancheMultiplication(
        alpha_n, alpha_p, W_cm);
    const double G_BTBT  = PhysicsEngine::kaneBTBT(
        ui.probeField, m.A_kane, m.B_kane, m.btbt_isDirect);

    ImGui::Text("alpha_n = %.3e cm^-1", alpha_n);
    ImGui::Text("alpha_p = %.3e cm^-1", alpha_p);
    if (std::isinf(Mfac)) {
        ImGui::TextColored(ImVec4(1.0f, 0.30f, 0.20f, 1.0f),
            "M -> infinity   (AVALANCHE BREAKDOWN)");
    } else {
        ImGui::Text("M       = %.3f", Mfac);
    }
    ImGui::Text("G_BTBT  = %.3e cm^-3 / s", G_BTBT);
    ImGui::TextDisabled("Kane: P=%s, A=%.2e, B=%.2e V/cm",
                        m.btbt_isDirect ? "2 (direct)" : "5/2 (indirect)",
                        m.A_kane, m.B_kane);

    // ---- Painter-grid coupling: re-evaluate at peak grid field ---------
    if (dd.deviceMode() == DeviceMode::Painter) {
        ImGui::SeparatorText("Grid evaluation (peak |E| of psi map)");
        const float E_peak = dd.peakElectricField();
        const double an = PhysicsEngine::chynowethRate(
            E_peak, m.alpha_inf_n, m.E_crit_n, m.chyn_m);
        const double ap = PhysicsEngine::chynowethRate(
            E_peak, m.alpha_inf_p, m.E_crit_p, m.chyn_m);
        const double Gb = PhysicsEngine::kaneBTBT(
            E_peak, m.A_kane, m.B_kane, m.btbt_isDirect);
        ImGui::Text("|E|_peak = %.3e V/cm", E_peak);
        ImGui::Text("alpha_n  = %.3e cm^-1", an);
        ImGui::Text("alpha_p  = %.3e cm^-1", ap);
        ImGui::Text("G_BTBT   = %.3e cm^-3 / s", Gb);
    }

    ImGui::End();
}


// =============================================================================
// 3D Topology window  [Phase 3]
// -----------------------------------------------------------------------------
// Industry-standard TCAD heatmap of the chosen scalar field over the
// device cross-section. Uses ImPlot::PlotHeatmap (always available). The
// underlying buffers are exposed by DriftDiffusion::psi() etc., so a
// future drop-in to ImPlot3D::PlotSurface only needs to swap the call
// site -- the data layout is already row-major (j*W + i).
//
// To minimise work per frame, derived fields (n, p, |E|) are built into
// a *static* scratch vector that grows once and is reused; this keeps
// the hot-path zero-allocation after the first call.
// =============================================================================
void drawTopology3DWindow(const DriftDiffusion& dd,
                          const PhysicsEngine& physics,
                          UIState& ui)
{
    if (!ImGui::Begin("3D Topology")) { ImGui::End(); return; }

    const char* fieldLabels[] = {
        "Electrostatic potential psi(x,y) [V]",
        "Electron density n(x,y) [cm^-3, log10]",
        "Hole density p(x,y) [cm^-3, log10]",
        "Electric field |E|(x,y) [V/cm, log10]",
        "Lattice Temperature T(x,y) [K]",          // [Phase 5]
        "Heat source H(x,y) [W/cm^3, log10]",      // [Phase 5]
        "Bandgap E_g(x,y,T) [eV]",                 // [Phase 5]
    };
    ImGui::Combo("Field", &ui.topo3DField,
                 fieldLabels, IM_ARRAYSIZE(fieldLabels));
    HelpMarker(
        "Industry-standard TCAD heatmap of a scalar field over the "
        "device cross-section. n/p/|E| are shown on log10 scale; psi is "
        "linear. Data is row-major (j*W+i) and is ready for an ImPlot3D::"
        "PlotSurface drop-in when that library is available.");

    const int W = dd.width();
    const int H = dd.height();
    if (W <= 0 || H <= 0) { ImGui::End(); return; }

    // Scratch buffer -- grown once, reused forever.
    static std::vector<float> scratch;
    scratch.resize(static_cast<std::size_t>(W) * static_cast<std::size_t>(H));

    const double n_i = physics.getIntrinsicCarrier();
    const double V_T = PhysicsEngine::thermalVoltage(physics.getTemperature());

    // Populate scratch according to the selection.
    auto safeLog10 = [](double v) {
        return std::log10(std::max(v, 1.0e-30));
    };
    switch (ui.topo3DField) {
        case 0: { // psi
            const auto& src = dd.psi();
            for (std::size_t k = 0; k < scratch.size(); ++k)
                scratch[k] = src[k];
            break;
        }
        case 1: { // n (log10)
            for (int j = 0; j < H; ++j)
                for (int i = 0; i < W; ++i)
                    scratch[static_cast<std::size_t>(j * W + i)] =
                        static_cast<float>(safeLog10(
                            dd.nDensityAt(i, j, n_i, V_T)));
            break;
        }
        case 2: { // p (log10)
            for (int j = 0; j < H; ++j)
                for (int i = 0; i < W; ++i)
                    scratch[static_cast<std::size_t>(j * W + i)] =
                        static_cast<float>(safeLog10(
                            dd.pDensityAt(i, j, n_i, V_T)));
            break;
        }
        case 3: { // |E| (log10)
            for (int j = 0; j < H; ++j)
                for (int i = 0; i < W; ++i) {
                    const double Em = dd.electricFieldMagAt(i, j);
                    scratch[static_cast<std::size_t>(j * W + i)] =
                        static_cast<float>(safeLog10(Em));
                }
            break;
        }
        case 4: { // [Phase 5] Lattice temperature T(x,y) [K] (linear)
            for (int j = 0; j < H; ++j)
                for (int i = 0; i < W; ++i)
                    scratch[static_cast<std::size_t>(j * W + i)] =
                        dd.temperatureFieldAt(i, j);
            break;
        }
        case 5: { // [Phase 5] Heat source H(x,y) [W/cm^3, log10]
            const auto& H_field = dd.heatGenField();
            for (int j = 0; j < H; ++j)
                for (int i = 0; i < W; ++i) {
                    const std::size_t k = static_cast<std::size_t>(j * W + i);
                    scratch[k] = static_cast<float>(
                        safeLog10(static_cast<double>(H_field[k])));
                }
            break;
        }
        case 6: { // [Phase 5] Local bandgap E_g(x,y,T) [eV] (linear)
            for (int j = 0; j < H; ++j)
                for (int i = 0; i < W; ++i)
                    scratch[static_cast<std::size_t>(j * W + i)] =
                        static_cast<float>(dd.localBandgapAt(i, j));
            break;
        }
        default: break;
    }

    // Auto-range.
    float vmin =  std::numeric_limits<float>::infinity();
    float vmax = -std::numeric_limits<float>::infinity();
    for (float v : scratch) {
        if (v < vmin) vmin = v;
        if (v > vmax) vmax = v;
    }
    if (vmin == vmax) { vmin -= 1.0f; vmax += 1.0f; }

    ImGui::Text("range: [%.4g, %.4g]", vmin, vmax);

    // ImPlot heatmap. PlotHeatmap takes (label, data, rows, cols, scale_min,
    // scale_max, fmt, bounds_min, bounds_max).
    if (ImPlot::BeginPlot("##topo_heatmap",
                          ImVec2(-1.0f, ImGui::GetContentRegionAvail().y - 24.f),
                          ImPlotFlags_NoLegend))
    {
        // Phase 5: hot fields use a perceptually-warm map so the eye
        // immediately reads them as "temperature".  Others stay viridis.
        const ImPlotColormap cmap = (ui.topo3DField == 4 || ui.topo3DField == 5)
            ? ImPlotColormap_Hot : ImPlotColormap_Viridis;
        ImPlot::PushColormap(cmap);
        ImPlot::SetupAxes("x [cell]", "y [cell]",
                          ImPlotAxisFlags_NoGridLines,
                          ImPlotAxisFlags_NoGridLines | ImPlotAxisFlags_Invert);
        ImPlot::SetupAxesLimits(0, W, 0, H);
        ImPlot::PlotHeatmap("field",
            scratch.data(), H, W,
            vmin, vmax,
            nullptr,
            ImPlotPoint(0, 0), ImPlotPoint(W, H));
        ImPlot::PopColormap();
        ImPlot::EndPlot();
    }

    ImGui::TextDisabled(
        "Future ImPlot3D upgrade: dd.psi().data(), rows=%d, cols=%d, "
        "row-major.", H, W);

    ImGui::End();
}


// =============================================================================
// Drift-Diffusion inspection
// =============================================================================
void drawCrystalInfoWindow(const DriftDiffusion& dd,
                           const PhysicsEngine& physics)
{
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

    // Bulk Shockley-Read-Hall recombination rate, computed from the global
    // n / p / n_i and the material's tau_n / tau_p:
    //
    //   R_SRH = (np - n_i^2) / [tau_p (n + n_i) + tau_n (p + n_i)]
    //
    // (Pierret Ch. 5, Sze Sec. 1.5.5). In thermal equilibrium np = n_i^2
    // and R_SRH = 0; under injection or illumination it becomes positive.
    {
        const auto& mat = physics.getMaterial();
        const double R = PhysicsEngine::recombSRH(
            physics.getTotalElectronConc(),
            physics.getTotalHoleConc(),
            physics.getIntrinsicCarrier(),
            mat.tau_n, mat.tau_p);
        ImGui::Text("R_SRH:      %.3e cm^-3 / s", R);
        HelpMarker(
            "Net Shockley-Read-Hall recombination rate using the material's "
            "tau_n, tau_p. Detailed balance: equals zero in thermal "
            "equilibrium, positive when np > n_i^2 (injection / illumination).");
    }

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

        // ---- Figures of merit derived from material + region geometry ----
        // Grid topology: emitter spans 22% of grid width; base 8%; the
        // remainder is collector. Physical width is taken as 1 um per cell
        // (standard textbook BJT scale).
        const auto& mat = physics.getMaterial();
        constexpr double cell_um = 1.0;
        const double W_E = dd.width() * 0.22 * cell_um * 1.0e-4;   // cm
        const double W_B = dd.width() * 0.08 * cell_um * 1.0e-4;
        constexpr double N_E = 1.0e19;     // typical heavy emitter
        constexpr double N_B = 1.0e16;     // typical base
        const double T = physics.getTemperature();
        const double D_n = mat.mu_L_n_300 * phys::k_B * T;   // Einstein
        const double D_p = mat.mu_L_p_300 * phys::k_B * T;
        const double L_n = std::sqrt(D_n * mat.tau_n);
        const double alpha_T_factor = std::clamp(
            1.0 - 0.5 * (W_B / L_n) * (W_B / L_n), 0.0, 0.99999);
        const double gamma_eff = PhysicsEngine::bjtEmitterEfficiency(
            D_n, N_E, W_E, D_p, N_B, W_B);
        const double beta_dc   = PhysicsEngine::bjtCurrentGain(
            gamma_eff, alpha_T_factor);
        const double early_fac = PhysicsEngine::earlyEffectFactor(
            dd.vCE(), mat.V_Early);

        ImGui::SeparatorText("BJT figures of merit");
        ImGui::Text("gamma:    %.5f",          gamma_eff);
        ImGui::Text("alpha_T:  %.5f",          alpha_T_factor);
        ImGui::Text("beta_DC:  %.1f",          beta_dc);
        ImGui::Text("V_Early:  %.1f V",        mat.V_Early);
        ImGui::Text("I_C/I_C0: %.3f (Early)",  early_fac);
        HelpMarker(
            "gamma = 1/(1 + D_p N_B W_E / D_n N_E W_B)   (emitter efficiency)\n"
            "alpha_T = 1 - W_B^2 / (2 L_n^2)             (transport factor)\n"
            "beta_DC = gamma alpha_T / (1 - gamma alpha_T)\n"
            "Early factor = 1 + V_CE / V_A");
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
    std::filesystem::path fontPath;
    (void)loadFont(font, fontPath);
    loadImGuiFont(fontPath);

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
        crystal->update(dt, *physics, *dd);
        dd->step(dt);
        physics->setDriftDiffusionExcess(dd->globalExcess());

        // ---- Self-consistent solver (Phase 1 Poisson or Phase 3 Gummel) -
        // Painter mode only -- procedural Bulk/BJT topologies have no
        // user-painted doping. Bias != 0 -> Gummel; otherwise plain
        // equilibrium Poisson. The user can also force Gummel ON via the
        // checkbox in Controls.
        if (dd->deviceMode() == DeviceMode::Painter) {
            const auto& mat = physics->getMaterial();
            const double n_i = physics->getIntrinsicCarrier();
            const double V_T =
                PhysicsEngine::thermalVoltage(physics->getTemperature());

            const bool useGummel = ui.gummelEnable
                || std::abs(ui.V_bias) > 1.0e-4f;

            if (useGummel) {
                if (ui.transientEnabled) {
                    // Transient (Backward Euler) -- AC probe & step
                    // pulse modify V_a inside stepTransient.
                    if (ui.stepPulseArmed && !ui.stepPulseFired) {
                        const float V_target = ui.V_bias + ui.stepPulseAmp;
                        dd->setAppliedBias(V_target);
                        ui.stepPulseFired = true;
                        ui.stepPulseT0    = dd->simulationTime();
                    }
                    ui.gummelResidual = dd->stepTransient(
                        n_i, V_T, mat.epsilon_r, mat,
                        physics->getTemperature(),
                        ui.gummelOuter,
                        ui.gummelPoissonInner,
                        ui.gummelContinuityInner,
                        ui.gummelOmegaPsi,
                        ui.gummelOmegaPhi);
                    if (ui.acProbeEnabled) ui.V_bias = dd->appliedBias();

                    // Phase 5: Wachutka local heat solve riding on the
                    // same dt as the carrier step. The result feeds back
                    // into n_i / E_g via refreshLocalParamCache() inside
                    // the next Gummel call -- closes the runaway loop.
                    if (ui.wachutkaThermal) {
                        ui.dT_step_last = dd->solveHeatEquation(
                            n_i, V_T,
                            physics->getElectronMobility(),
                            physics->getHoleMobility(),
                            dd->timeStep());
                        ui.T_peak_local = dd->maxTemperature();
                        ui.T_mean_local = dd->meanTemperature();
                    }
                } else {
                    ui.gummelResidual = dd->solveGummel(
                        n_i, V_T, mat.epsilon_r,
                        physics->getElectronMobility(),
                        physics->getHoleMobility(),
                        mat,
                        physics->getTemperature(),
                        ui.gummelOuter,
                        ui.gummelPoissonInner,
                        ui.gummelContinuityInner,
                        ui.gummelOmegaPsi,
                        ui.gummelOmegaPhi);
                }
                ui.poissonResidual = ui.gummelResidual;
                ui.lastJ_terminal  = dd->terminalCurrentDensity(
                    n_i, V_T,
                    physics->getElectronMobility(),
                    physics->getHoleMobility());

                // Append to oscilloscope ring (cap at ~3000 samples).
                ui.oscTime.push_back(
                    static_cast<float>(dd->simulationTime()));
                ui.oscCurrent.push_back(
                    static_cast<float>(ui.lastJ_terminal));
                if (ui.oscTime.size() > 3000) {
                    ui.oscTime.erase   (ui.oscTime.begin());
                    ui.oscCurrent.erase(ui.oscCurrent.begin());
                }
            } else if (ui.poissonAutoSolve) {
                ui.poissonResidual = dd->solvePoisson(
                    n_i, V_T, mat.epsilon_r,
                    ui.poissonIters,
                    ui.poissonOmega);
            }
        }

        // Re-roll dopants/carriers periodically when the doping state could
        // have changed (cheap; only resamples carrier visualization).
        // Done implicitly by Controls window callbacks via rebuild path.
        // (Here we just refresh carrier population every frame so optical
        //  and drift-diffusion changes are visible.)
        if (static int frame = 0; (++frame % 30) == 0)
            crystal->rebuild(*physics);

        // ---- Render off-screen targets ------------------------------------
        crystal->render(*physics, *dd, ui.heatmapMode, ui.showVectorField);

        // BandView: spatial bent bands when in Painter mode (or when the
        // user explicitly opts in); fall back to the legacy flat diagram
        // for Bulk / BJT modes where there is no spatially varying psi.
        bands->setMode(ui.bandSpatialMode ? BandView::Mode::Spatial
                                          : BandView::Mode::Flat);
        if (bands->mode() == BandView::Mode::Spatial) {
            bands->renderSpatial(*physics, *dd, ui.bandCutRow, font);
        } else {
            bands->render(*physics, font);
        }

        // ---- ImGui frame --------------------------------------------------
        ImGui::SFML::Update(window, delta);

        (void)beginDockspaceHost(requestResetLayout);
        drawMenuBar(running, ui, *physics, *dd, requestResetLayout);
        ImGui::End();   // close DockHost

        drawControlsWindow     (*physics, *dd, ui);
        drawReadoutsWindow     (*physics, *dd);
        drawLiveOscilloscope   (*physics, *dd, ui, dt);
        drawSpectrumWindow     (*physics);
        drawIVCurveWindow      (*dd);
        drawCrystalViewWindow  (*crystal, *dd, *physics, ui);
        drawDevicePainterWindow(*dd, *physics, ui);
        drawBandViewWindow     (*bands, *dd, ui);
        drawBreakdownWindow    (*physics, *dd, ui);
        drawTopology3DWindow   (*dd, *physics, ui);
        drawCrystalInfoWindow  (*dd, *physics);

        // ---- Final composite ---------------------------------------------
        window.clear(palette::WindowBg);
        ImGui::SFML::Render(window);
        window.display();
    }

    ImPlot::DestroyContext();
    ImGui::SFML::Shutdown();
    return 0;
}
