#pragma once

// =============================================================================
// DriftDiffusion.hpp -- Phase 7: interactive 2D carrier diffusion
// -----------------------------------------------------------------------------
// Solves the minority-carrier continuity equation on a coarse 2D grid:
//
//       d n / d t  =  D nabla^2 n  +  G(x, y, t)  -  n / tau
//
// using an explicit FTCS (Forward-Time Centered-Space) finite-difference
// scheme. The user can click on the crystal view to deposit a Gaussian
// "laser spot" (or "heat point") source which then diffuses outward and
// recombines on the lifetime time scale tau.
//
// Pedagogical, not quantitative: D and tau are scaled so that one click is
// visible for several seconds at 60 fps without any unit-system gymnastics.
//
// References: Sze Ch. 1, Pierret "Semiconductor Device Fundamentals" Ch. 6.
// =============================================================================

#include <cstddef>
#include <vector>


class DriftDiffusion {
public:
    // Grid resolution chosen to comfortably fit inside the CrystalView's
    // RenderTexture without making the FTCS step too expensive (60x40 cells
    // is ~2400 grid points, evaluated each frame in well under 1 ms).
    DriftDiffusion(int gridW = 60, int gridH = 40);

    // Material lifetime / diffusion-coefficient scale (called when material
    // changes; we re-tune the integrator constants per material so e.g. GaAs
    // diffuses faster than Si visually).
    void configureForMaterial(double D_scale, double tau_seconds);

    // ---- Sources --------------------------------------------------------
    // Drop a Gaussian source (intensity, sigma) at normalized image coords
    // (u, v) in [0, 1]. Sources persist until clear() is called.
    void addSource(float u, float v, float intensity = 1.0f, float sigma = 0.04f);

    // Drop everything.
    void clear();

    // ---- Time stepping --------------------------------------------------
    void step(float dt_seconds);

    // ---- Inspection -----------------------------------------------------
    [[nodiscard]] int   width()      const noexcept { return m_W; }
    [[nodiscard]] int   height()     const noexcept { return m_H; }
    [[nodiscard]] float at(int i, int j) const noexcept;     // density [a.u.]
    [[nodiscard]] float maxValue()   const noexcept;
    [[nodiscard]] float meanValue()  const noexcept;

    // Read-only access to the entire grid for renderers.
    [[nodiscard]] const std::vector<float>& data() const noexcept { return m_n; }

    // ---- Approximate "average excess" returned to the PhysicsEngine ----
    // Used so the global readouts (sigma, R_H, ...) reflect the existence
    // of the user-placed sources. Scaled to a physical-looking dN [cm^-3].
    [[nodiscard]] double globalExcess() const noexcept;

private:
    // Grid bookkeeping
    int  m_W;
    int  m_H;
    std::vector<float> m_n;       // current carrier excess
    std::vector<float> m_n_next;  // double-buffer
    std::vector<float> m_G;       // persistent source distribution

    // Integrator coefficients (visualisation units, not physical cm^2/s).
    float m_D    = 0.18f;   // diffusion strength per step
    float m_tau  = 1.6f;    // recombination time constant [s]
    float m_dN_scale = 1.0e15f; // converts unit-less density -> cm^-3 readout

    // Helpers
    [[nodiscard]] inline std::size_t idx(int i, int j) const noexcept {
        return static_cast<std::size_t>(j) * static_cast<std::size_t>(m_W)
             + static_cast<std::size_t>(i);
    }
};
