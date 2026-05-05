#pragma once

// =============================================================================
// DriftDiffusion.hpp -- Phase 6: 2D coupled electrothermal solver
// -----------------------------------------------------------------------------
// What started in Phase 7 of the SFML era as a single carrier-diffusion grid
// now hosts three coupled fields plus a region map:
//
//      +-------------------------------------------------------------+
//      |  n(x,y)  : excess minority-carrier density   [a.u.]         |
//      |  T(x,y)  : lattice temperature                [K]           |
//      |  G(x,y)  : persistent generation source       [a.u./s]      |
//      |  region  : Bulk / Emitter / Base / Collector  (BJT mode)    |
//      +-------------------------------------------------------------+
//
// Two PDEs are stepped each frame using FTCS with a per-step CFL clamp:
//
//   d n/d t  =  D_n nabla^2 n + G(x,y) - n / tau                 (carriers)
//   d T/d t  =  alpha_T nabla^2 T + H_gen / (rho Cp)             (heat)
//
//   where  alpha_T = kappa / (rho Cp)   is the thermal diffusivity and
//   H_gen has two contributions:
//       Joule:          H_J  ~  sigma * E^2     (V_CE^2 / L^2 in BJT)
//       Recombination:  H_R  ~  R_SRH * E_g     (excess n / tau)
//
// Time-scale mismatch is handled by sub-stepping: M carrier substeps and K
// thermal substeps per render frame, with M >= K because carriers diffuse
// faster (D_n >> alpha_T).
//
// All buffers are pre-allocated in the constructor; per-frame execution is
// strictly zero-allocation.
// =============================================================================

#include <cstddef>
#include <cstdint>
#include <vector>

#include "Material.hpp"


// -----------------------------------------------------------------------------
// Top-level device mode. Bulk = a single homogeneous wafer with optional
// laser sources. NpnBjt = N-type emitter + thin P-base + N-collector with
// V_BE / V_CE biasing (textbook NPN).
// -----------------------------------------------------------------------------
enum class DeviceMode : std::uint8_t { Bulk = 0, NpnBjt = 1 };


// -----------------------------------------------------------------------------
// Per-cell region tag. Bulk cells participate in carrier diffusion only;
// emitter/collector cells additionally enforce Dirichlet boundary conditions
// for n every step (forward injection / sweep-out).
// -----------------------------------------------------------------------------
enum class CellRegion : std::uint8_t {
    Bulk      = 0,
    Emitter   = 1,
    Base      = 2,
    Collector = 3,
};


class DriftDiffusion {
public:
    DriftDiffusion(int gridW = 60, int gridH = 40);

    // ---- Material / configuration ----------------------------------------
    // Re-tunes the integrator constants (D_n, alpha_T) for the new material.
    // Should be called whenever PhysicsEngine::setMaterial is called.
    void configureForMaterial(const material::Profile& mat);

    // ---- Source placement (laser / heat point) ---------------------------
    void addSource(float u, float v,
                   float intensity = 1.0f,
                   float sigma     = 0.04f);
    void clear();

    // ---- Thermal --------------------------------------------------------
    void setAmbientTemperature(float T_kelvin);
    [[nodiscard]] float ambientTemperature() const noexcept { return m_T_ambient; }
    [[nodiscard]] float temperatureAt(int i, int j) const noexcept;
    [[nodiscard]] float maxTemperature()  const noexcept;
    [[nodiscard]] float meanTemperature() const noexcept;

    // ---- Device mode (Bulk vs BJT) ---------------------------------------
    void setDeviceMode(DeviceMode mode);
    [[nodiscard]] DeviceMode deviceMode() const noexcept { return m_mode; }
    void setBjtVoltages(float V_BE, float V_CE);
    [[nodiscard]] float vBE() const noexcept { return m_VBE; }
    [[nodiscard]] float vCE() const noexcept { return m_VCE; }
    [[nodiscard]] CellRegion regionAt(int i, int j) const noexcept;

    // Approximate collector-emitter current (sweep-out rate). Useful for
    // reporting transistor "operating point" in the readouts panel.
    [[nodiscard]] float collectorCurrent() const noexcept { return m_I_C; }

    // ---- Time stepping --------------------------------------------------
    // Advances both fields by dt_seconds (real wall time). Internally:
    //   * stepCarriers is sub-stepped M=4 times.
    //   * stepThermal  is sub-stepped K=1 time.
    void step(float dt_seconds);

    // ---- Inspection -----------------------------------------------------
    [[nodiscard]] int   width()      const noexcept { return m_W; }
    [[nodiscard]] int   height()     const noexcept { return m_H; }
    [[nodiscard]] float at(int i, int j) const noexcept;       // n grid
    [[nodiscard]] float maxValue()   const noexcept;
    [[nodiscard]] float meanValue()  const noexcept;
    [[nodiscard]] const std::vector<float>& data()       const noexcept { return m_n; }
    [[nodiscard]] const std::vector<float>& temperature()const noexcept { return m_T; }
    [[nodiscard]] const std::vector<std::uint8_t>& regions() const noexcept { return m_region; }

    // Average excess carrier density translated into a "cm^-3" readout for
    // the global PhysicsEngine UI.
    [[nodiscard]] double globalExcess() const noexcept;

    // Average temperature change above ambient (for thermal readouts).
    [[nodiscard]] float deltaTaverage() const noexcept;

private:
    // ---- Substep kernels ------------------------------------------------
    void stepCarriers(float dt);
    void stepThermal (float dt);
    void applyBjtBoundaries();
    void rebuildRegionMap();

    // ---- Index helpers --------------------------------------------------
    [[nodiscard]] inline std::size_t idx(int i, int j) const noexcept {
        return static_cast<std::size_t>(j) * static_cast<std::size_t>(m_W)
             + static_cast<std::size_t>(i);
    }

    // ---- Geometry -------------------------------------------------------
    int m_W;
    int m_H;

    // ---- Carrier field (double-buffered) --------------------------------
    std::vector<float> m_n;
    std::vector<float> m_n_next;

    // ---- Source field ---------------------------------------------------
    std::vector<float> m_G;

    // ---- Thermal field --------------------------------------------------
    std::vector<float> m_T;
    std::vector<float> m_T_next;

    // ---- Region map (BJT) -----------------------------------------------
    std::vector<std::uint8_t> m_region;

    // ---- Material-tunable integrator constants --------------------------
    float m_D_n        = 0.04f;   // carrier diffusion strength per unit dt
    float m_alpha_T    = 0.018f;  // thermal diffusivity (rescaled)
    float m_tau        = 1.6f;    // carrier recombination lifetime [s]
    float m_dN_scale   = 1.0e15f; // visualisation -> cm^-3 scaling
    float m_T_ambient  = 300.0f;  // [K]

    // ---- Device-mode state ---------------------------------------------
    DeviceMode m_mode = DeviceMode::Bulk;
    float m_VBE = 0.0f;     // base-emitter forward bias [V]
    float m_VCE = 0.0f;     // collector-emitter bias    [V]
    float m_I_C = 0.0f;     // last-frame collector current (a.u.)

    // Cached material thermal source coupling
    const material::Profile* m_material = nullptr;
};
