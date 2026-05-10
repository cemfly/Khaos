#pragma once

// =============================================================================
// DriftDiffusion.hpp -- 2D coupled electrothermal solver
//
//   Author : dex / cemfly-april2026
//   License: MIT
// -----------------------------------------------------------------------------
// Three coupled scalar fields plus a region map share one grid:
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
#include <span>
#include <vector>

#include "Material.hpp"


// -----------------------------------------------------------------------------
// Top-level device mode. Bulk = a single homogeneous wafer with optional
// laser sources. NpnBjt = N-type emitter + thin P-base + N-collector with
// V_BE / V_CE biasing (textbook NPN). Painter = user-defined topology (the
// painter UI stamps Nd/Na onto the grid; rebuildRegionMap leaves it alone).
// -----------------------------------------------------------------------------
enum class DeviceMode : std::uint8_t { Bulk = 0, NpnBjt = 1, Painter = 2 };


// -----------------------------------------------------------------------------
// Per-cell region tag. Bulk cells participate in carrier diffusion only;
// emitter/collector cells additionally enforce Dirichlet boundary conditions
// for n every step. NDoped/PDoped are Painter-mode tags applied by the
// brush; the Poisson solver does not care about the tag (it reads Nd, Na
// directly), but the visualisation does.
// -----------------------------------------------------------------------------
enum class CellRegion : std::uint8_t {
    Bulk      = 0,
    Emitter   = 1,
    Base      = 2,
    Collector = 3,
    NDoped    = 4,
    PDoped    = 5,
};


// -----------------------------------------------------------------------------
// Brush selection used by the Device Painter UI.
// -----------------------------------------------------------------------------
enum class BrushKind : std::uint8_t {
    None    = 0,
    NDopant = 1,    // stamps donors    (Nd += dose)
    PDopant = 2,    // stamps acceptors (Na += dose)
    Eraser  = 3,    // resets the painted cell to intrinsic
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

    // ---- Device painter [Phase 1] ---------------------------------------
    //
    // Stamps a circular brush of radius `radius_cells` around UV (0..1)
    // at concentration `dose` [cm^-3]. Strictly zero-allocation; iterates
    // a clipped square AABB. NDopant/PDopant are additive; Eraser zeroes.
    void paintBrush(float u, float v,
                    BrushKind kind,
                    double dose      = 1.0e17,
                    int radius_cells = 3) noexcept;
    void clearDoping() noexcept;

    [[nodiscard]] double netDopingAt(int i, int j) const noexcept;
    [[nodiscard]] double donorAt    (int i, int j) const noexcept;
    [[nodiscard]] double acceptorAt (int i, int j) const noexcept;
    void setCellPitchCm(float h) noexcept;
    [[nodiscard]] float cellPitchCm() const noexcept { return m_cell_pitch_cm; }

    // ---- Poisson solver  [Phase 1] --------------------------------------
    //
    // Solves equilibrium  eps_s grad^2 psi = -q (p(psi) - n(psi) + Nd - Na)
    // by Gauss-Seidel under-relaxation. Returns L2 residual of the last
    // sweep (in volts). All buffers preallocated; zero-allocation hot-path.
    [[nodiscard]] double solvePoisson(double n_i_cm3,
                                      double V_T,
                                      double epsilon_r,
                                      int    iterations = 80,
                                      double omega      = 0.85) noexcept;

    [[nodiscard]] float psiAt(int i, int j) const noexcept;
    [[nodiscard]] float bandShiftAt(int i, int j) const noexcept;

    // Electric field magnitude at (i,j) from psi gradient.
    //   E = -grad(psi);   |E| = sqrt(E_x^2 + E_y^2)         [V/cm]
    // Used by Chynoweth / Kane BTBT plotting.
    [[nodiscard]] float electricFieldMagAt(int i, int j) const noexcept;
    [[nodiscard]] float electricFieldX  (int i, int j) const noexcept;
    [[nodiscard]] float electricFieldY  (int i, int j) const noexcept;
    [[nodiscard]] float peakElectricField() const noexcept;

    // ---- Spatial cut sampling [Phase 2: BandView] -----------------------
    //
    // Fills caller-owned buffers along a horizontal cut at j = j_row:
    //   psi[i], Ec[i], Ev[i]                 -- length = m_W
    //
    // Pre-allocated callers; we never resize. Returns false if the row is
    // out of range or the spans are too small.
    [[nodiscard]] bool sampleHorizontalCut(
        int j_row,
        std::span<float> out_psi,
        std::span<float> out_Ec,
        std::span<float> out_Ev,
        double Ec0, double Ev0) const noexcept;

    // Sample electric-field magnitude along a horizontal cut.
    [[nodiscard]] bool sampleEFieldCut(
        int j_row, std::span<float> out_E) const noexcept;

    // Sample net doping along a horizontal cut (signed Nd-Na).
    [[nodiscard]] bool sampleDopingCut(
        int j_row, std::span<float> out_net) const noexcept;

    // ---- Phase 3 -- Non-equilibrium / Gummel solver ---------------------
    //
    // External bias applied between the painted contacts (anode at +V_a,
    // cathode at 0).  Detected automatically: leftmost / rightmost columns
    // are scanned for majority dopant; the side with more acceptors is
    // tagged anode, the side with more donors is tagged cathode. The user
    // does not need to mark contacts explicitly.
    void  setAppliedBias(float V_a) noexcept;
    [[nodiscard]] float appliedBias() const noexcept { return m_V_bias; }

    // Quasi-Fermi state. In equilibrium both buffers are zero everywhere,
    // and the equations reduce to the Boltzmann / Phase-1 limit.
    [[nodiscard]] float phiN(int i, int j) const noexcept;
    [[nodiscard]] float phiP(int i, int j) const noexcept;

    // Local density evaluated from psi, phi_n / phi_p:
    //   n = n_i exp((psi - phi_n) / V_T),
    //   p = n_i exp((phi_p - psi) / V_T).
    [[nodiscard]] double nDensityAt(int i, int j,
                                    double n_i, double V_T) const noexcept;
    [[nodiscard]] double pDensityAt(int i, int j,
                                    double n_i, double V_T) const noexcept;

    // Gummel iteration: outer loop Poisson -> Continuity_n -> Continuity_p.
    // Each inner solver is Gauss-Seidel under-relaxation on a 5-point
    // stencil. mu_n / mu_p are the *bulk* low-field mobilities supplied
    // by the engine (Matthiessen / Arora).  Returns the L2 change in psi
    // across the last outer iteration (volts) -- a useful convergence
    // indicator for the UI.
    //
    // Strictly zero-allocation: every working buffer is preallocated at
    // construction.  Scharfetter-Gummel-style flux discretization is
    // approximated here with arithmetic-mean face densities; for the
    // pedagogical real-time UI this is stable up to ~10 V reverse bias
    // and ~1 V forward bias on a 1 um grid.
    //
    // Reference: Gummel, IEEE TED-11 (1964) 455; Selberherr Ch. 6;
    // Sze Sec. 2.3.
    //
    // T_lattice [K] is required to evaluate the Matthiessen mobility at
    // each face (lattice + impurity scattering depend on T).
    [[nodiscard]] double solveGummel(
        double n_i, double V_T, double epsilon_r,
        double mu_n_bulk, double mu_p_bulk,
        const material::Profile& mat,
        double T_lattice         = 300.0,
        int    outer_iters       = 6,
        int    poisson_inner     = 30,
        int    continuity_inner  = 20,
        double omega_psi         = 0.85,
        double omega_phi         = 1.00) noexcept;

    // ---- Phase 3 -- Quasi-Fermi cuts (BandView spatial render) ---------
    [[nodiscard]] bool sampleQuasiFermiCut(
        int j_row,
        std::span<float> out_phi_n,
        std::span<float> out_phi_p) const noexcept;

    // Average current density across the device (mid-column proxy):
    // J = J_n + J_p averaged along i = W/2.  Pierret App. C.
    [[nodiscard]] double terminalCurrentDensity(
        double n_i, double V_T,
        double mu_n, double mu_p) const noexcept;

    // ---- Bulk data accessors (3D heatmap window, future ImPlot3D) ------
    [[nodiscard]] const std::vector<float>& psi()         const noexcept { return m_psi; }
    [[nodiscard]] const std::vector<float>& phiNField()   const noexcept { return m_phi_n; }
    [[nodiscard]] const std::vector<float>& phiPField()   const noexcept { return m_phi_p; }
    [[nodiscard]] const std::vector<float>& donorField()  const noexcept { return m_Nd; }
    [[nodiscard]] const std::vector<float>& acceptorField()const noexcept{ return m_Na; }
    [[nodiscard]] const std::vector<float>& nDensField()  const noexcept { return m_n_dens; }
    [[nodiscard]] const std::vector<float>& pDensField()  const noexcept { return m_p_dens; }

    // ---- Phase 4 -- Transient (Backward Euler) -------------------------
    //
    // Time-stepping API. dt > 0 enables the transient term in continuity:
    //
    //   (n_new - n_old) / dt - div(J_n)/q = G - R
    //
    // Backward Euler is unconditionally stable -- no CFL constraint on dt.
    // Pass 0 (or call setTransientEnabled(false)) to fall back to pure
    // steady-state Gummel.
    void  setTimeStep(double dt_s) noexcept;
    [[nodiscard]] double timeStep() const noexcept { return m_dt; }
    void  setTransientEnabled(bool e) noexcept { m_transient_enabled = e; }
    [[nodiscard]] bool transientEnabled() const noexcept { return m_transient_enabled; }
    [[nodiscard]] double simulationTime() const noexcept { return m_sim_time; }
    void  resetSimulationTime() noexcept { m_sim_time = 0.0; }

    // Advances the simulation by m_dt: copies n,p -> n_old,p_old, applies
    // any AC probe to V_a, then runs one Gummel iteration. Returns the
    // L2 psi residual (same convention as solveGummel).
    [[nodiscard]] double stepTransient(
        double n_i, double V_T, double epsilon_r,
        const material::Profile& mat,
        double T_lattice,
        int    outer_iters       = 4,
        int    poisson_inner     = 20,
        int    continuity_inner  = 12,
        double omega_psi         = 0.85,
        double omega_phi         = 1.00) noexcept;

    // ---- Phase 4 -- AC probe (small-signal stimulus) -------------------
    //
    //   V_a(t) = V_DC + amp * sin(2 pi f t)
    //
    // The DC base is the value passed to setAppliedBias before AC is
    // turned on; the AC envelope rides on top. Useful with stepTransient
    // to visualise capacitive charging on Live Oscilloscope.
    void  setACProbe(bool enabled, double freq_Hz, double amp_V) noexcept;
    [[nodiscard]] bool   acEnabled() const noexcept { return m_ac_enabled; }
    [[nodiscard]] double acFreq()    const noexcept { return m_ac_freq; }
    [[nodiscard]] double acAmp()     const noexcept { return m_ac_amp; }

    // ---- Phase 4 -- Small-signal estimators (DC perturbation) ----------
    //
    //  G_DC = dJ/dV  : run two extra Gummel solves at V +/- dV, finite-
    //                  difference the terminal current.  Returns S/cm^2.
    //  C_dep / C_diff : reported via PhysicsEngine helpers using a
    //                  detected depletion width and tau.
    [[nodiscard]] double smallSignalConductance(
        double n_i, double V_T, double epsilon_r,
        const material::Profile& mat,
        double T_lattice,
        double dV = 0.005) noexcept;

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

    // ---- Phase 3 inner solvers ------------------------------------------
    //
    // The shared core of Phase-1 solvePoisson and Phase-3 solveGummel:
    // performs `sweeps` Gauss-Seidel iterations of the nonlinear Poisson
    // equation reading m_phi_n, m_phi_p (zero -> equilibrium reduction).
    [[nodiscard]] double solvePoissonInner(double n_i, double V_T,
                                           double epsilon_r,
                                           int sweeps,
                                           double omega) noexcept;

    // Scharfetter-Gummel continuity solvers [Phase 4]
    //
    //   J_n(i+1/2) = (q V_T / h) mu_face [B(x) n_{i+1} - B(-x) n_i]
    //   J_p(i+1/2) = -(q V_T / h) mu_face [B(-x) p_{i+1} - B(x) p_i]
    //   x = (psi_{i+1} - psi_i) / V_T
    //
    // Per-face mu evaluated from local doping (Matthiessen) and field
    // (Caughey-Thomas).  Backward-Euler time term added when m_dt > 0
    // and m_transient_enabled.  Stable up to >50 V reverse bias.
    //
    // mu_n_bulk / mu_p_bulk are the engine's bulk fallbacks (used at
    // boundaries where doping is intrinsic and Matthiessen would return
    // very high values).
    //
    // Reference: Selberherr Sec. 6.2 / Scharfetter-Gummel IEEE-ED 16 (1969).
    void solveContinuityElectron(double n_i, double V_T,
                                 double mu_n_bulk,
                                 const material::Profile& mat,
                                 double T_lattice,
                                 int sweeps, double omega) noexcept;
    void solveContinuityHole    (double n_i, double V_T,
                                 double mu_p_bulk,
                                 const material::Profile& mat,
                                 double T_lattice,
                                 int sweeps, double omega) noexcept;

    // Apply Dirichlet BCs at the outer columns (anode / cathode), set
    // by majority-dopant detection. Called at the start of each Gummel
    // outer iteration.
    void applyContactBoundaries(double n_i, double V_T) noexcept;

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

    // ---- Region map (BJT / Painter) -------------------------------------
    std::vector<std::uint8_t> m_region;

    // ---- Doping fields (painter-driven) [Phase 1] -----------------------
    std::vector<float> m_Nd;          // donor    concentration [cm^-3]
    std::vector<float> m_Na;          // acceptor concentration [cm^-3]

    // ---- Electrostatics [Phase 1] ---------------------------------------
    std::vector<float> m_psi;         // electrostatic potential [V]
    std::vector<float> m_psi_next;    // double buffer for sweeps

    // Cell pitch in cm. Default 1 um is a textbook visualisation scale;
    // for sharper junctions reduce to ~L_D / 4 (Si N=1e16 -> ~10 nm).
    float m_cell_pitch_cm = 1.0e-4f;

    // ---- Phase 3 -- Non-equilibrium quasi-Fermi state -------------------
    std::vector<float> m_phi_n;       // electron quasi-Fermi potential [V]
    std::vector<float> m_phi_p;       // hole     quasi-Fermi potential [V]
    float              m_V_bias = 0.0f;   // applied anode bias [V]

    // Cached contact roles (computed once per Gummel call, not per sweep).
    // 0 = none, 1 = anode, 2 = cathode. Stored as a per-cell uint8 so we
    // can tag arbitrary user-painted layouts (not just the outer columns).
    std::vector<std::uint8_t> m_contact;
    bool                      m_contacts_dirty = true;

    // ---- Phase 4 -- Density buffers (SG continuity primary unknowns) ---
    std::vector<float> m_n_dens;      // electron density [cm^-3]
    std::vector<float> m_p_dens;      // hole     density [cm^-3]
    std::vector<float> m_n_dens_old;  // previous time step (BE transient)
    std::vector<float> m_p_dens_old;
    bool               m_dens_seeded   = false;   // initial Boltzmann seed?

    // ---- Phase 4 -- Transient state ------------------------------------
    double m_dt                = 0.0;     // time step [s]; 0 = steady-state
    double m_sim_time          = 0.0;     // accumulated time [s]
    bool   m_transient_enabled = false;

    // ---- Phase 4 -- AC small-signal probe ------------------------------
    bool   m_ac_enabled = false;
    double m_ac_freq    = 1.0e6;          // [Hz]
    double m_ac_amp     = 0.005;          // [V]
    float  m_V_dc_base  = 0.0f;           // DC operating point used by AC

    // Internal: refresh m_phi_n / m_phi_p from m_n_dens / m_p_dens / m_psi.
    // Called once per Gummel iteration so the BandView quasi-Fermi cuts
    // stay consistent with the SG-iterated densities.
    void refreshQuasiFermiFromDensities(double n_i, double V_T) noexcept;

    // Internal: bootstrap n,p from local Boltzmann limit using current
    // psi / phi buffers. Cheap O(N); first solveGummel call seeds, later
    // calls reuse the converged state.
    void seedDensitiesFromBoltzmann(double n_i, double V_T) noexcept;

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
