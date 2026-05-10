// =============================================================================
// Preset.cpp -- versioned JSON save/load for the simulator state
//
//   Author : dex / cemfly-april2026
//   License: MIT
// -----------------------------------------------------------------------------
// One-file dump of every user-controlled knob in Khaos: physics inputs,
// painter doping map, device mode, bias, transient + AC settings.
// Layout documented in Preset.hpp.
//
// Round-trip strategy
//   * Save  : pull every value via existing const-accessors, write JSON.
//   * Load  : parse fully *first* (so a half-bad file leaves state alone),
//             then apply settings in a single pass.
// =============================================================================

#include "Preset.hpp"

#include <chrono>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>

#include <nlohmann/json.hpp>

#include "DriftDiffusion.hpp"
#include "Material.hpp"
#include "PhysicsEngine.hpp"


using json = nlohmann::json;


namespace preset {

// =============================================================================
// Error -> string (UI status banner)
// =============================================================================
std::string toString(Error e) noexcept {
    switch (e) {
        case Error::None:              return "ok";
        case Error::OpenFailed:        return "cannot open file";
        case Error::ParseFailed:       return "invalid JSON";
        case Error::SchemaMismatch:    return "unsupported preset schema";
        case Error::FieldMissing:      return "required field missing";
        case Error::FieldTypeMismatch: return "wrong field type";
        case Error::SizeMismatch:      return "doping array size mismatch";
        case Error::WriteFailed:       return "write failed";
    }
    return "unknown error";
}


// =============================================================================
// Local helpers
// =============================================================================
namespace {

// Enum <-> string converters. Centralised so the same vocabulary appears
// on both save and load and a typo crashes the load (rather than silently
// picking a default).
[[nodiscard]] const char* materialKindToStr(material::Kind k) noexcept {
    switch (k) {
        case material::Kind::Silicon:    return "Si";
        case material::Kind::GaAs:       return "GaAs";
        case material::Kind::Germanium:  return "Ge";
    }
    return "Si";
}
[[nodiscard]] std::expected<material::Kind, Error>
materialKindFromStr(std::string_view s) noexcept {
    if (s == "Si")   return material::Kind::Silicon;
    if (s == "GaAs") return material::Kind::GaAs;
    if (s == "Ge")   return material::Kind::Germanium;
    return std::unexpected(Error::FieldTypeMismatch);
}

[[nodiscard]] const char* dopingTypeToStr(DopingType t) noexcept {
    switch (t) {
        case DopingType::Intrinsic: return "Intrinsic";
        case DopingType::NType:     return "NType";
        case DopingType::PType:     return "PType";
    }
    return "Intrinsic";
}
[[nodiscard]] std::expected<DopingType, Error>
dopingTypeFromStr(std::string_view s) noexcept {
    if (s == "Intrinsic") return DopingType::Intrinsic;
    if (s == "NType")     return DopingType::NType;
    if (s == "PType")     return DopingType::PType;
    return std::unexpected(Error::FieldTypeMismatch);
}

[[nodiscard]] const char* mobilityModelToStr(MobilityModel m) noexcept {
    return (m == MobilityModel::Matthiessen) ? "Matthiessen" : "Arora";
}
[[nodiscard]] std::expected<MobilityModel, Error>
mobilityModelFromStr(std::string_view s) noexcept {
    if (s == "Matthiessen") return MobilityModel::Matthiessen;
    if (s == "Arora")       return MobilityModel::Arora;
    return std::unexpected(Error::FieldTypeMismatch);
}

[[nodiscard]] const char* deviceModeToStr(DeviceMode m) noexcept {
    switch (m) {
        case DeviceMode::Bulk:    return "Bulk";
        case DeviceMode::NpnBjt:  return "NpnBjt";
        case DeviceMode::Painter: return "Painter";
    }
    return "Bulk";
}
[[nodiscard]] std::expected<DeviceMode, Error>
deviceModeFromStr(std::string_view s) noexcept {
    if (s == "Bulk")    return DeviceMode::Bulk;
    if (s == "NpnBjt")  return DeviceMode::NpnBjt;
    if (s == "Painter") return DeviceMode::Painter;
    return std::unexpected(Error::FieldTypeMismatch);
}

// UTC ISO-8601 timestamp for the "captured_at" field. Lightweight; no
// strftime needed in this <chrono>-only C++20 path.
[[nodiscard]] std::string utcTimestamp() {
    const auto now = std::chrono::system_clock::now();
    const auto t   = std::chrono::system_clock::to_time_t(now);
    std::tm tm{};
#if defined(_WIN32)
    gmtime_s(&tm, &t);
#else
    gmtime_r(&t, &tm);
#endif
    std::ostringstream oss;
    oss << (tm.tm_year + 1900) << '-'
        << std::setw(2) << std::setfill('0') << (tm.tm_mon + 1) << '-'
        << std::setw(2) << std::setfill('0') << tm.tm_mday      << 'T'
        << std::setw(2) << std::setfill('0') << tm.tm_hour      << ':'
        << std::setw(2) << std::setfill('0') << tm.tm_min       << ':'
        << std::setw(2) << std::setfill('0') << tm.tm_sec       << 'Z';
    return oss.str();
}

} // namespace


// =============================================================================
// Save
// =============================================================================
std::expected<void, Error> save(
    const std::filesystem::path& path,
    const PhysicsEngine&         physics,
    const DriftDiffusion&        dd)
{
    json j;
    j["schema"]      = kSchemaMajor;
    j["schema_minor"]= kSchemaMinor;
    j["tool"]        = "Khaos / dex / cemfly-april2026";
    j["captured_at"] = utcTimestamp();

    // ---- Physics ---------------------------------------------------------
    json& jp = j["physics"];
    jp["material"]              = materialKindToStr(physics.getMaterialKind());
    jp["temperature_K"]         = physics.getTemperature();
    jp["doping_type"]           = dopingTypeToStr(physics.getDopingType());
    jp["doping_concentration"]  = physics.getDopingConcentration();
    jp["incomplete_ionization"] = physics.getIncompleteIonization();
    jp["mobility_model"]        = mobilityModelToStr(physics.getMobilityModel());
    jp["wavelength_nm"]         = physics.getWavelengthNm();
    jp["optical_enabled"]       = physics.getOpticalEnabled();
    jp["magnetic_field_T"]      = physics.getMagneticField();

    // ---- Device ----------------------------------------------------------
    json& jd = j["device"];
    jd["mode"]          = deviceModeToStr(dd.deviceMode());
    jd["V_BE"]          = dd.vBE();
    jd["V_CE"]          = dd.vCE();
    jd["V_bias"]        = dd.appliedBias();
    jd["ambient_T_K"]   = dd.ambientTemperature();
    jd["cell_pitch_cm"] = dd.cellPitchCm();
    jd["grid_w"]        = dd.width();
    jd["grid_h"]        = dd.height();

    // Painter doping arrays -- store as compact JSON arrays of double.
    // Float would suffice for size, but JSON has no distinct float type;
    // double is honest and round-trips cleanly through serialisation.
    // nlohmann/json's vector adapter takes ownership and serialises
    // directly -- avoids the json::array_t internal-type assumption.
    const auto& Nd = dd.donorField();
    const auto& Na = dd.acceptorField();
    std::vector<double> NdD; NdD.reserve(Nd.size());
    std::vector<double> NaD; NaD.reserve(Na.size());
    for (float v : Nd) NdD.push_back(static_cast<double>(v));
    for (float v : Na) NaD.push_back(static_cast<double>(v));
    jd["Nd"] = std::move(NdD);
    jd["Na"] = std::move(NaD);

    // ---- Transient + AC --------------------------------------------------
    json& jt = j["transient"];
    jt["enabled"]    = dd.transientEnabled();
    jt["dt_s"]       = dd.timeStep();
    jt["sim_time_s"] = dd.simulationTime();

    json jac;
    jac["enabled"] = dd.acEnabled();
    jac["freq_Hz"] = dd.acFreq();
    jac["amp_V"]   = dd.acAmp();
    jt["ac"]       = std::move(jac);

    // ---- Write -----------------------------------------------------------
    std::ofstream out(path, std::ios::binary);
    if (!out) return std::unexpected(Error::OpenFailed);
    out << j.dump(2);       // 2-space indent: human-diffable; small overhead
    if (!out.good()) return std::unexpected(Error::WriteFailed);
    return {};
}


// =============================================================================
// Load
// =============================================================================
namespace {

// Tiny field-extraction helpers that map JSON errors onto our Error
// enum and short-circuit propagation via std::expected. Avoids the
// nlohmann::json::exception clutter at call sites.
template <typename T>
[[nodiscard]] std::expected<T, Error> required(const json& j,
                                               const char* key) noexcept
{
    auto it = j.find(key);
    if (it == j.end()) return std::unexpected(Error::FieldMissing);
    try { return it->get<T>(); }
    catch (...) { return std::unexpected(Error::FieldTypeMismatch); }
}

} // namespace

std::expected<std::string, Error> load(
    const std::filesystem::path& path,
    PhysicsEngine&               physics,
    DriftDiffusion&              dd)
{
    std::ifstream in(path, std::ios::binary);
    if (!in) return std::unexpected(Error::OpenFailed);

    json j;
    try { in >> j; }
    catch (const json::parse_error&) {
        return std::unexpected(Error::ParseFailed);
    }

    // ---- Schema check ----------------------------------------------------
    auto schema = required<int>(j, "schema");
    if (!schema) return std::unexpected(schema.error());
    if (*schema != kSchemaMajor)
        return std::unexpected(Error::SchemaMismatch);

    auto jp = j.find("physics");
    auto jd = j.find("device");
    if (jp == j.end() || jd == j.end())
        return std::unexpected(Error::FieldMissing);

    // ---- Parse all fields first (no mutation yet) ------------------------
    auto mat_s   = required<std::string>(*jp, "material");
    auto T       = required<double>     (*jp, "temperature_K");
    auto dt_s    = required<std::string>(*jp, "doping_type");
    auto N       = required<double>     (*jp, "doping_concentration");
    auto inc     = required<bool>       (*jp, "incomplete_ionization");
    auto mob_s   = required<std::string>(*jp, "mobility_model");
    auto lam     = required<double>     (*jp, "wavelength_nm");
    auto opt     = required<bool>       (*jp, "optical_enabled");
    auto B_T     = required<double>     (*jp, "magnetic_field_T");

    auto dev_s   = required<std::string>(*jd, "mode");
    auto Vbe     = required<double>     (*jd, "V_BE");
    auto Vce     = required<double>     (*jd, "V_CE");
    auto Vbias   = required<double>     (*jd, "V_bias");
    auto Tamb    = required<double>     (*jd, "ambient_T_K");
    auto pitch   = required<double>     (*jd, "cell_pitch_cm");
    auto gw      = required<int>        (*jd, "grid_w");
    auto gh      = required<int>        (*jd, "grid_h");

    if (!mat_s || !T || !dt_s || !N || !inc || !mob_s || !lam || !opt || !B_T
        || !dev_s || !Vbe || !Vce || !Vbias || !Tamb || !pitch || !gw || !gh)
        return std::unexpected(Error::FieldMissing);

    auto kind     = materialKindFromStr(*mat_s);
    auto doping   = dopingTypeFromStr  (*dt_s);
    auto mobility = mobilityModelFromStr(*mob_s);
    auto devmode  = deviceModeFromStr  (*dev_s);
    if (!kind || !doping || !mobility || !devmode)
        return std::unexpected(Error::FieldTypeMismatch);

    // Grid size compatibility: we don't resize the live DriftDiffusion
    // (its buffers are pre-allocated for performance), so the saved grid
    // *must* match the runtime grid.
    if (*gw != dd.width() || *gh != dd.height())
        return std::unexpected(Error::SizeMismatch);

    // Doping arrays
    auto itNd = jd->find("Nd");
    auto itNa = jd->find("Na");
    if (itNd == jd->end() || itNa == jd->end())
        return std::unexpected(Error::FieldMissing);
    if (!itNd->is_array() || !itNa->is_array())
        return std::unexpected(Error::FieldTypeMismatch);

    const std::size_t expected_cells =
        static_cast<std::size_t>(*gw) * static_cast<std::size_t>(*gh);
    if (itNd->size() != expected_cells || itNa->size() != expected_cells)
        return std::unexpected(Error::SizeMismatch);

    // Transient + AC (optional in v1; absent = default off)
    bool   transEnabled = false;
    double dt_s_t       = 1.0e-9;
    bool   acEnabled    = false;
    double acFreq       = 1.0e6;
    double acAmp        = 0.005;
    if (auto jt = j.find("transient"); jt != j.end()) {
        if (auto x = required<bool>  (*jt, "enabled")) transEnabled = *x;
        if (auto x = required<double>(*jt, "dt_s"))    dt_s_t       = *x;
        if (auto jac = jt->find("ac"); jac != jt->end()) {
            if (auto x = required<bool>  (*jac, "enabled")) acEnabled = *x;
            if (auto x = required<double>(*jac, "freq_Hz")) acFreq    = *x;
            if (auto x = required<double>(*jac, "amp_V"))   acAmp     = *x;
        }
    }

    // ---- Apply (everything past this point should not fail) --------------
    physics.setMaterial            (*kind);
    physics.setTemperature         (*T);
    physics.setDopingType          (*doping);
    physics.setDopingConcentration (*N);
    physics.setIncompleteIonization(*inc);
    physics.setMobilityModel       (*mobility);
    physics.setWavelengthNm        (*lam);
    physics.setOpticalEnabled      (*opt);
    physics.setMagneticField       (*B_T);

    dd.setAmbientTemperature(static_cast<float>(*Tamb));
    dd.setCellPitchCm       (static_cast<float>(*pitch));

    // Repaint the device map. clearDoping zeroes both Nd/Na and resets
    // region tags; setDopingAt then rewrites each cell.
    dd.clearDoping();
    for (int j2 = 0; j2 < *gh; ++j2) {
        for (int i = 0; i < *gw; ++i) {
            const std::size_t k = static_cast<std::size_t>(j2 * (*gw) + i);
            dd.setDopingAt(i, j2,
                (*itNd)[k].get<double>(),
                (*itNa)[k].get<double>());
        }
    }

    dd.setDeviceMode    (*devmode);
    dd.setBjtVoltages   (static_cast<float>(*Vbe),
                         static_cast<float>(*Vce));
    dd.setAppliedBias   (static_cast<float>(*Vbias));
    dd.setTimeStep      (dt_s_t);
    dd.setTransientEnabled(transEnabled);
    dd.setACProbe       (acEnabled, acFreq, acAmp);

    std::ostringstream summary;
    summary << "preset v" << *schema << ".x loaded (grid "
            << *gw << "x" << *gh << ", "
            << materialKindToStr(*kind) << ", V_a = "
            << *Vbias << " V)";
    return summary.str();
}

} // namespace preset
