#pragma once

// =============================================================================
// Preset.hpp -- versioned JSON save/load for the simulator state
//
//   Author : dex / cemfly-april2026
//   License: MIT
// -----------------------------------------------------------------------------
// Serialises the complete user-visible simulation state to a single
// `.khaos.json` file using the nlohmann/json schema below. Designed
// around three principles:
//
//   1. Self-contained -- one file, no external references.
//   2. Versioned     -- a top-level "schema" field lets future loaders
//                       migrate or reject older / unknown captures.
//   3. Lossless on round-trip for everything the user sees in the UI.
//
// Schema (v1):
//   {
//     "schema": 1,
//     "tool":   "Khaos / dex / cemfly-april2026",
//     "captured_at": "<ISO-8601 UTC>",
//     "physics": {
//        "material":            "Si" | "GaAs" | "Ge",
//        "temperature_K":       300.0,
//        "doping_type":         "Intrinsic" | "NType" | "PType",
//        "doping_concentration": 1.0e16,
//        "incomplete_ionization": false,
//        "mobility_model":      "Matthiessen" | "Arora",
//        "wavelength_nm":       1500.0,
//        "optical_enabled":     false,
//        "magnetic_field_T":    0.0
//     },
//     "device": {
//        "mode":         "Bulk" | "NpnBjt" | "Painter",
//        "V_BE":         0.0,
//        "V_CE":         0.0,
//        "V_bias":       0.0,
//        "ambient_T_K":  300.0,
//        "cell_pitch_cm": 1.0e-4,
//        "grid_w":       60,
//        "grid_h":       40,
//        "Nd":           [...flat row-major, length = grid_w*grid_h...],
//        "Na":           [...]
//     },
//     "transient": {
//        "enabled": false,
//        "dt_s":    1.0e-9,
//        "sim_time_s": 0.0,
//        "ac": { "enabled": false, "freq_Hz": 1.0e6, "amp_V": 0.005 }
//     }
//   }
//
// Note: the live psi / n / p / T fields are intentionally NOT saved -- the
// solver re-derives them by running Gummel on the captured Nd / Na / V_a
// map at load time. That keeps preset files small (just the topology)
// and prevents stale solutions from masking schema drift.
// =============================================================================

#include <expected>
#include <filesystem>
#include <string>


class PhysicsEngine;
class DriftDiffusion;


namespace preset {

// -----------------------------------------------------------------------------
// Schema constants -- bump major (incompatible) or minor (additive) on
// every change. Loaders compare against this and reject anything they
// cannot interpret.
// -----------------------------------------------------------------------------
inline constexpr int kSchemaMajor = 1;
inline constexpr int kSchemaMinor = 0;


// -----------------------------------------------------------------------------
// Error reporting -- a tiny enum + free-form detail string suit a UI
// status banner better than exceptions. std::expected lets the caller
// branch on success without try/catch overhead.
// -----------------------------------------------------------------------------
enum class Error : std::uint8_t {
    None              = 0,
    OpenFailed        = 1,    // file could not be opened (permission, path)
    ParseFailed       = 2,    // not valid JSON
    SchemaMismatch    = 3,    // schema field missing / unsupported major
    FieldMissing      = 4,    // required field absent
    FieldTypeMismatch = 5,    // wrong JSON type
    SizeMismatch      = 6,    // doping array size != grid_w * grid_h
    WriteFailed       = 7,    // file could not be written
};

[[nodiscard]] std::string toString(Error e) noexcept;

struct LoadResult {
    bool        loaded = false;
    std::string detail;
};


// -----------------------------------------------------------------------------
// Save:  serialise (physics + dd) to a UTF-8 JSON file at `path`.
// Returns std::expected: error code carries the failure reason.
// -----------------------------------------------------------------------------
[[nodiscard]] std::expected<void, Error> save(
    const std::filesystem::path& path,
    const PhysicsEngine&         physics,
    const DriftDiffusion&        dd);

// -----------------------------------------------------------------------------
// Load:  read a JSON file and apply it to (physics, dd). The two objects
// are mutated in-place; if loading fails partway, the caller's state is
// left untouched (we parse fully first, then apply).
// Returns a human-readable detail line on success (e.g. "preset v1.0 OK").
// -----------------------------------------------------------------------------
[[nodiscard]] std::expected<std::string, Error> load(
    const std::filesystem::path& path,
    PhysicsEngine&               physics,
    DriftDiffusion&              dd);

} // namespace preset
