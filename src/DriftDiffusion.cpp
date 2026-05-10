// =============================================================================
// DriftDiffusion.cpp
//
//   Author : dex / cemfly-april2026
//   License: MIT
// =============================================================================

#include "DriftDiffusion.hpp"

#include <algorithm>
#include <cmath>

#include "PhysicsEngine.hpp"   // SRH / Auger / Kane statics for continuity src


// =============================================================================
// Construction -- pre-allocate every grid buffer once and never again.
// =============================================================================
DriftDiffusion::DriftDiffusion(int gridW, int gridH)
    : m_W(gridW), m_H(gridH)
    , m_n       (static_cast<std::size_t>(gridW * gridH), 0.0f)
    , m_n_next  (static_cast<std::size_t>(gridW * gridH), 0.0f)
    , m_G       (static_cast<std::size_t>(gridW * gridH), 0.0f)
    , m_T       (static_cast<std::size_t>(gridW * gridH), 300.0f)
    , m_T_next  (static_cast<std::size_t>(gridW * gridH), 300.0f)
    , m_region  (static_cast<std::size_t>(gridW * gridH),
                 static_cast<std::uint8_t>(CellRegion::Bulk))
    , m_Nd      (static_cast<std::size_t>(gridW * gridH), 0.0f)
    , m_Na      (static_cast<std::size_t>(gridW * gridH), 0.0f)
    , m_psi     (static_cast<std::size_t>(gridW * gridH), 0.0f)
    , m_psi_next(static_cast<std::size_t>(gridW * gridH), 0.0f)
    , m_phi_n   (static_cast<std::size_t>(gridW * gridH), 0.0f)
    , m_phi_p   (static_cast<std::size_t>(gridW * gridH), 0.0f)
    , m_contact (static_cast<std::size_t>(gridW * gridH), 0u)
    // Phase 4: density buffers (SG continuity primary unknowns)
    , m_n_dens     (static_cast<std::size_t>(gridW * gridH), 0.0f)
    , m_p_dens     (static_cast<std::size_t>(gridW * gridH), 0.0f)
    , m_n_dens_old (static_cast<std::size_t>(gridW * gridH), 0.0f)
    , m_p_dens_old (static_cast<std::size_t>(gridW * gridH), 0.0f)
{}


// =============================================================================
// Material / integrator tuning
//
//   D_n       : carrier diffusion strength (Si baseline = 0.04 dimensionless).
//               Scales linearly with the lattice mobility ratio so GaAs
//               diffuses ~6x faster than Si in the visualisation.
//   alpha_T   : kappa / (rho Cp). Si gets ~0.018 per-frame coefficient at
//               default tuning; other materials scale relative to this.
// =============================================================================
void DriftDiffusion::configureForMaterial(const material::Profile& mat) {
    m_material = &mat;

    const double D_scale = mat.mu_L_n_300 / 1414.0;  // Si baseline
    m_D_n     = static_cast<float>(0.04 * D_scale);

    // Thermal diffusivity in physical units: alpha = kappa / rho_cp [cm^2/s].
    // We squash it to a per-frame coefficient by dividing by a pedagogical
    // length scale^2; the Si value lands in a comfortable visual range
    // (~0.018) and the other materials track relatively.
    const double alpha_phys  = mat.kappa / mat.rho_cp;        // cm^2/s
    constexpr double alpha_si = 1.50 / 1.65;                  // ~0.91
    m_alpha_T = static_cast<float>(0.018 * (alpha_phys / alpha_si));

    m_tau = 1.6f;  // visualisation lifetime
}


// =============================================================================
// Source helpers
// =============================================================================
void DriftDiffusion::addSource(float u, float v, float intensity, float sigma) {
    u = std::clamp(u, 0.0f, 1.0f);
    v = std::clamp(v, 0.0f, 1.0f);
    const float i0  = u * (m_W - 1);
    const float j0  = v * (m_H - 1);
    const float two_sigma_sq_inv = 1.0f / (2.0f * sigma * sigma * (m_W * m_H));

    for (int j = 0; j < m_H; ++j) {
        for (int i = 0; i < m_W; ++i) {
            const float du = (i - i0);
            const float dv = (j - j0);
            const float r2 = du * du + dv * dv;
            const float w  = std::exp(-r2 * two_sigma_sq_inv);
            m_G[idx(i, j)] += intensity * w;
        }
    }
}

void DriftDiffusion::clear() {
    std::fill(m_n.begin(),      m_n.end(),      0.0f);
    std::fill(m_n_next.begin(), m_n_next.end(), 0.0f);
    std::fill(m_G.begin(),      m_G.end(),      0.0f);
    std::fill(m_T.begin(),      m_T.end(),      m_T_ambient);
    std::fill(m_T_next.begin(), m_T_next.end(), m_T_ambient);
}


// =============================================================================
// Ambient temperature -- written to grid edges every thermal step (Dirichlet).
// =============================================================================
void DriftDiffusion::setAmbientTemperature(float T) {
    m_T_ambient = std::clamp(T, 50.0f, 1200.0f);
}


// =============================================================================
// Device mode
// =============================================================================
void DriftDiffusion::setDeviceMode(DeviceMode m) {
    if (m_mode == m) return;
    m_mode = m;
    // Painter mode owns the region map -- never overwrite user paint.
    if (m == DeviceMode::Bulk || m == DeviceMode::NpnBjt) {
        rebuildRegionMap();
    }
}

void DriftDiffusion::setBjtVoltages(float V_BE, float V_CE) {
    m_VBE = std::clamp(V_BE, 0.0f, 1.2f);
    m_VCE = std::clamp(V_CE, 0.0f, 5.0f);
}

CellRegion DriftDiffusion::regionAt(int i, int j) const noexcept {
    if (i < 0 || j < 0 || i >= m_W || j >= m_H) return CellRegion::Bulk;
    return static_cast<CellRegion>(m_region[idx(i, j)]);
}


// -----------------------------------------------------------------------------
// Build the BJT region map: |  E  | B |     C     |
//
//   * Emitter   : leftmost ~22% of width  (heavily doped, electron source)
//   * Base      : narrow strip ~8%         (P-type, sweep-through region)
//   * Collector : remaining ~70%           (electron sink)
//
// Bulk mode wipes everything back to CellRegion::Bulk.
// -----------------------------------------------------------------------------
void DriftDiffusion::rebuildRegionMap() {
    if (m_mode == DeviceMode::Painter) return;   // user-owned topology
    if (m_mode == DeviceMode::Bulk) {
        std::fill(m_region.begin(), m_region.end(),
                  static_cast<std::uint8_t>(CellRegion::Bulk));
        return;
    }

    const int x_emitter_end = static_cast<int>(m_W * 0.22f);
    const int x_base_end    = static_cast<int>(m_W * 0.30f);

    for (int j = 0; j < m_H; ++j) {
        for (int i = 0; i < m_W; ++i) {
            CellRegion r;
            if      (i <  x_emitter_end) r = CellRegion::Emitter;
            else if (i <  x_base_end)    r = CellRegion::Base;
            else                          r = CellRegion::Collector;
            m_region[idx(i, j)] = static_cast<std::uint8_t>(r);
        }
    }
}


// =============================================================================
// Carrier substep -- FTCS with optional BJT Dirichlet boundary post-step
//
//   n_next[i,j] = n[i,j]
//               + D_eff * (n[i+1,j] + n[i-1,j] + n[i,j+1] + n[i,j-1] - 4 n[i,j])
//               + G[i,j] * dt
//               - n[i,j] * dt / tau
//
//   D_eff is clamped to <= 0.24 per FTCS step regardless of dt.
//   Reflective Neumann at walls except where the BJT contacts override it.
// =============================================================================
void DriftDiffusion::stepCarriers(float dt) {
    const float D_eff = std::min(0.24f, m_D_n * dt * 60.0f);
    const float decay = dt / m_tau;

    for (int j = 0; j < m_H; ++j) {
        for (int i = 0; i < m_W; ++i) {
            const int ip = std::min(i + 1, m_W - 1);
            const int im = std::max(i - 1, 0);
            const int jp = std::min(j + 1, m_H - 1);
            const int jm = std::max(j - 1, 0);

            const float c   = m_n[idx(i,  j )];
            const float lap = m_n[idx(ip, j )] + m_n[idx(im, j )]
                            + m_n[idx(i,  jp)] + m_n[idx(i,  jm)]
                            - 4.0f * c;

            float v = c + D_eff * lap + m_G[idx(i, j)] * dt - c * decay;
            v = std::max(0.0f, v);
            m_n_next[idx(i, j)] = v;
        }
    }
    m_n.swap(m_n_next);

    if (m_mode == DeviceMode::NpnBjt) applyBjtBoundaries();
}


// =============================================================================
// BJT boundary conditions (Dirichlet)
//
//   Emitter   : n = N_inj * exp(q V_BE / kT)        (forward injection)
//   Collector : n = 0                               (reverse-biased sink)
//
// Also accumulates the collector-emitter "current" as the sum of the change
// across the collector boundary -- a proxy for sweep-out rate.
// =============================================================================
void DriftDiffusion::applyBjtBoundaries() {
    // kT at the *current* ambient temperature -- so the V_BE response is
    // properly temperature-dependent (cold devices need more V_BE for the
    // same injection level; hot devices saturate faster).
    const float kT_eV      = 8.617333262e-5f * m_T_ambient;
    const float thermalArg = std::clamp(m_VBE / kT_eV, 0.0f, 25.0f);
    const float n_emitter  = 0.6f * std::exp(thermalArg);
    // n_emitter capped because the visualisation grid is unitless; full
    // exp(40) would saturate immediately and hide structure. The clamp
    // exposes the steep V_BE response without going off-screen.

    // Collector sink strength scales mildly with V_CE (more reverse bias
    // -> more aggressive depletion) and provides a useful BJT current proxy.
    const float sweep_factor = 0.0f + m_VCE * 0.03f;

    float swept = 0.0f;
    for (int j = 0; j < m_H; ++j) {
        for (int i = 0; i < m_W; ++i) {
            const auto r = static_cast<CellRegion>(m_region[idx(i, j)]);
            if (r == CellRegion::Emitter) {
                m_n[idx(i, j)] = n_emitter;
            } else if (r == CellRegion::Collector) {
                const float prev = m_n[idx(i, j)];
                swept += prev * (1.0f - sweep_factor);
                m_n[idx(i, j)] *= sweep_factor;     // sink (close to 0)
            }
        }
    }
    m_I_C = swept;   // proportional to collector current
}


// =============================================================================
// Thermal substep -- FTCS for the heat equation with Dirichlet edges
//
//   T_next[i,j] = T[i,j]
//               + alpha_eff * (Laplace stencil)
//               + H_gen[i,j] / (rho Cp) * dt
//
//   Edges are clamped to ambient T every step (heat sink contacts) so the
//   peak temperature is bounded.
//
//   H_gen has two pieces:
//     * Joule        : Joule heating from local current. In bulk mode it
//                      tracks the local carrier excess (laser injection
//                      heats the spot); in BJT mode it ramps with V_CE^2
//                      along the active path between emitter and collector.
//     * Recombination: minority-carrier recombination dumps E_g of energy
//                      per pair, scaled by lifetime.
// =============================================================================
void DriftDiffusion::stepThermal(float dt) {
    const float alpha_eff = std::min(0.24f, m_alpha_T * dt * 60.0f);
    const float Eg_J  = (m_material ? static_cast<float>(m_material->Eg0) : 1.12f);
    const float rhoCp = (m_material ? static_cast<float>(m_material->rho_cp) : 1.65f);

    // Joule pre-factor for BJT: more current -> more dissipation.
    // For bulk we let the source field G act as a "laser absorption" heater.
    const bool  bjt = (m_mode == DeviceMode::NpnBjt);
    const float joule_bjt = bjt
        ? std::min(40.0f, m_VCE * m_VCE * 4.0f) : 0.0f;

    for (int j = 0; j < m_H; ++j) {
        for (int i = 0; i < m_W; ++i) {
            // Dirichlet edges (always reset to ambient).
            if (i == 0 || j == 0 || i == m_W - 1 || j == m_H - 1) {
                m_T_next[idx(i, j)] = m_T_ambient;
                continue;
            }

            const float c   = m_T[idx(i,  j )];
            const float lap = m_T[idx(i+1, j  )] + m_T[idx(i-1, j  )]
                            + m_T[idx(i,   j+1)] + m_T[idx(i,   j-1)]
                            - 4.0f * c;

            // Heat generation per unit volume per unit time (a.u.)
            float H_gen = 0.0f;

            // -- Recombination heat -----------------------------------------
            //    R_SRH ~ n / tau ; each pair releases ~E_g of energy.
            const float n_local = m_n[idx(i, j)];
            H_gen += (n_local / m_tau) * Eg_J * 0.4f;

            // -- Joule heat (mode-specific) --------------------------------
            if (bjt) {
                const auto r = static_cast<CellRegion>(m_region[idx(i, j)]);
                if (r == CellRegion::Base || r == CellRegion::Collector)
                    H_gen += n_local * joule_bjt;
            } else {
                // Bulk: laser absorption heats the spot.
                H_gen += m_G[idx(i, j)] * 6.0f;
            }

            float T_new = c + alpha_eff * lap + (H_gen / rhoCp) * dt;

            // Clamp to a sane physical envelope so a runaway never breaks
            // the visualisation -- but allow plenty of headroom for users
            // to *see* the runaway happen first.
            T_new = std::clamp(T_new, 50.0f, 1200.0f);
            m_T_next[idx(i, j)] = T_new;
        }
    }
    m_T.swap(m_T_next);
}


// =============================================================================
// Master step -- multi-rate operator splitting
//
//   Carriers: M sub-steps  (fast diffusion, short tau)
//   Thermal : K sub-steps  (slow diffusion, long thermal time-constant)
// =============================================================================
void DriftDiffusion::step(float dt_seconds) {
    constexpr int M_carrier = 4;
    constexpr int K_thermal = 1;

    if (dt_seconds <= 0.0f) return;

    const float dt_n = dt_seconds / static_cast<float>(M_carrier);
    for (int s = 0; s < M_carrier; ++s) stepCarriers(dt_n);

    const float dt_T = dt_seconds / static_cast<float>(K_thermal);
    for (int s = 0; s < K_thermal; ++s) stepThermal(dt_T);
}


// =============================================================================
// Inspection helpers
// =============================================================================
float DriftDiffusion::at(int i, int j) const noexcept {
    if (i < 0 || j < 0 || i >= m_W || j >= m_H) return 0.0f;
    return m_n[idx(i, j)];
}

float DriftDiffusion::temperatureAt(int i, int j) const noexcept {
    if (i < 0 || j < 0 || i >= m_W || j >= m_H) return m_T_ambient;
    return m_T[idx(i, j)];
}

float DriftDiffusion::maxValue() const noexcept {
    return m_n.empty() ? 0.0f : *std::max_element(m_n.begin(), m_n.end());
}

float DriftDiffusion::meanValue() const noexcept {
    if (m_n.empty()) return 0.0f;
    double s = 0.0;
    for (float v : m_n) s += v;
    return static_cast<float>(s / m_n.size());
}

float DriftDiffusion::maxTemperature() const noexcept {
    return m_T.empty() ? m_T_ambient : *std::max_element(m_T.begin(), m_T.end());
}

float DriftDiffusion::meanTemperature() const noexcept {
    if (m_T.empty()) return m_T_ambient;
    double s = 0.0;
    for (float v : m_T) s += v;
    return static_cast<float>(s / m_T.size());
}

double DriftDiffusion::globalExcess() const noexcept {
    return static_cast<double>(meanValue()) * static_cast<double>(m_dN_scale);
}

float DriftDiffusion::deltaTaverage() const noexcept {
    return meanTemperature() - m_T_ambient;
}


// =============================================================================
// Painter [Phase 1]
// -----------------------------------------------------------------------------
// Brush is a circular stamp; iteration is bounded to the AABB of the brush
// radius so cost is O(radius^2). Strictly zero-allocation -- only touches
// pre-allocated vectors.
// =============================================================================
void DriftDiffusion::paintBrush(float u, float v,
                                BrushKind kind,
                                double    dose,
                                int       radius_cells) noexcept
{
    if (kind == BrushKind::None) return;
    u = std::clamp(u, 0.0f, 1.0f);
    v = std::clamp(v, 0.0f, 1.0f);

    const int ic = static_cast<int>(u * (m_W - 1));
    const int jc = static_cast<int>(v * (m_H - 1));
    const int r2 = radius_cells * radius_cells;

    const int i0 = std::max(0,        ic - radius_cells);
    const int i1 = std::min(m_W - 1,  ic + radius_cells);
    const int j0 = std::max(0,        jc - radius_cells);
    const int j1 = std::min(m_H - 1,  jc + radius_cells);

    const float dose_f = static_cast<float>(dose);

    for (int j = j0; j <= j1; ++j) {
        for (int i = i0; i <= i1; ++i) {
            const int di = i - ic;
            const int dj = j - jc;
            if (di * di + dj * dj > r2) continue;

            const std::size_t k = idx(i, j);
            switch (kind) {
                case BrushKind::NDopant:
                    m_Nd[k]    += dose_f;
                    m_region[k] = static_cast<std::uint8_t>(CellRegion::NDoped);
                    break;
                case BrushKind::PDopant:
                    m_Na[k]    += dose_f;
                    m_region[k] = static_cast<std::uint8_t>(CellRegion::PDoped);
                    break;
                case BrushKind::Eraser:
                    m_Nd[k]    = 0.0f;
                    m_Na[k]    = 0.0f;
                    m_region[k] = static_cast<std::uint8_t>(CellRegion::Bulk);
                    break;
                default: break;
            }
        }
    }
    m_contacts_dirty = true;
}

void DriftDiffusion::clearDoping() noexcept {
    std::fill(m_Nd.begin(),         m_Nd.end(),         0.0f);
    std::fill(m_Na.begin(),         m_Na.end(),         0.0f);
    std::fill(m_psi.begin(),        m_psi.end(),        0.0f);
    std::fill(m_psi_next.begin(),   m_psi_next.end(),   0.0f);
    std::fill(m_phi_n.begin(),      m_phi_n.end(),      0.0f);
    std::fill(m_phi_p.begin(),      m_phi_p.end(),      0.0f);
    std::fill(m_n_dens.begin(),     m_n_dens.end(),     0.0f);
    std::fill(m_p_dens.begin(),     m_p_dens.end(),     0.0f);
    std::fill(m_n_dens_old.begin(), m_n_dens_old.end(), 0.0f);
    std::fill(m_p_dens_old.begin(), m_p_dens_old.end(), 0.0f);
    std::fill(m_contact.begin(),    m_contact.end(),    0u);
    if (m_mode == DeviceMode::Painter) {
        std::fill(m_region.begin(), m_region.end(),
                  static_cast<std::uint8_t>(CellRegion::Bulk));
    }
    m_contacts_dirty = true;
    m_dens_seeded    = false;
    m_sim_time       = 0.0;
}

double DriftDiffusion::netDopingAt(int i, int j) const noexcept {
    if (i < 0 || j < 0 || i >= m_W || j >= m_H) return 0.0;
    const auto k = idx(i, j);
    return static_cast<double>(m_Nd[k]) - static_cast<double>(m_Na[k]);
}
double DriftDiffusion::donorAt   (int i, int j) const noexcept {
    if (i < 0 || j < 0 || i >= m_W || j >= m_H) return 0.0;
    return static_cast<double>(m_Nd[idx(i, j)]);
}
double DriftDiffusion::acceptorAt(int i, int j) const noexcept {
    if (i < 0 || j < 0 || i >= m_W || j >= m_H) return 0.0;
    return static_cast<double>(m_Na[idx(i, j)]);
}

void DriftDiffusion::setCellPitchCm(float h) noexcept {
    m_cell_pitch_cm = std::clamp(h, 1.0e-7f, 1.0e-2f);
}


// =============================================================================
// Self-consistent Poisson (equilibrium, Boltzmann statistics) [Phase 1]
// -----------------------------------------------------------------------------
// 5-point Laplacian on uniform pitch h:
//
//   psi_C^new = (1/4) * [ psi_E + psi_W + psi_N + psi_S
//                       + (q h^2 / eps_s) (p - n + Nd - Na) ]
//
// Under-relaxation (omega < 1) damps the strong nonlinearity from the
// exponential dependence of n,p on psi. Boundary: zero-derivative Neumann
// (ghost-cell mirroring) so the painted region drifts to its own
// equilibrium without spurious image charges.
// References: Selberherr Ch. 5; Sze Sec. 2.2.
// =============================================================================
double DriftDiffusion::solvePoisson(double n_i_cm3,
                                    double V_T,
                                    double epsilon_r,
                                    int    iterations,
                                    double omega) noexcept
{
    // Public Phase-1 entry point. Reads m_phi_n / m_phi_p (0 in pure
    // equilibrium -> identical to the original Boltzmann formulation).
    // Refresh contact BCs first so stale Gummel state can't leak in
    // (Phase 4 added contact tagging that solvePoissonInner now obeys).
    applyContactBoundaries(n_i_cm3, V_T);
    return solvePoissonInner(n_i_cm3, V_T, epsilon_r, iterations, omega);
}

double DriftDiffusion::solvePoissonInner(double n_i_cm3,
                                         double V_T,
                                         double epsilon_r,
                                         int    sweeps,
                                         double omega) noexcept
{
    constexpr double eps_0_local = 8.854187817e-14;     // F/cm
    const double eps_s    = eps_0_local * epsilon_r;
    constexpr double q    = 1.602176634e-19;
    const double h        = static_cast<double>(m_cell_pitch_cm);
    const double rho_coef = q * h * h / eps_s;          // [V * cm^3]

    if (V_T <= 0.0) return 0.0;

    double residual_l2 = 0.0;

    for (int it = 0; it < sweeps; ++it) {
        residual_l2 = 0.0;

        for (int j = 0; j < m_H; ++j) {
            for (int i = 0; i < m_W; ++i) {
                const std::size_t k = idx(i, j);

                // Skip cells fixed by Dirichlet contacts.
                if (m_contact[k] != 0u) continue;

                const int ip = std::min(i + 1, m_W - 1);
                const int im = std::max(i - 1, 0);
                const int jp = std::min(j + 1, m_H - 1);
                const int jm = std::max(j - 1, 0);

                const double psi_C = m_psi[k];
                const double psi_E = m_psi[idx(ip, j)];
                const double psi_W = m_psi[idx(im, j)];
                const double psi_N = m_psi[idx(i, jp)];
                const double psi_S = m_psi[idx(i, jm)];

                // Quasi-Fermi-aware densities (Phase 3). In equilibrium
                // m_phi_n = m_phi_p = 0 -> reduces to the Boltzmann form.
                const double phin_C = m_phi_n[k];
                const double phip_C = m_phi_p[k];
                const double xn = std::clamp((psi_C - phin_C) / V_T, -40.0, 40.0);
                const double xp = std::clamp((phip_C - psi_C) / V_T, -40.0, 40.0);
                const double n_C = n_i_cm3 * std::exp(xn);
                const double p_C = n_i_cm3 * std::exp(xp);

                const double rho =
                    p_C - n_C
                    + static_cast<double>(m_Nd[k])
                    - static_cast<double>(m_Na[k]);

                const double psi_new = 0.25 * (
                    psi_E + psi_W + psi_N + psi_S + rho_coef * rho
                );

                const double psi_relaxed = psi_C + omega * (psi_new - psi_C);
                const double dpsi        = psi_relaxed - psi_C;
                residual_l2             += dpsi * dpsi;

                m_psi[k] = static_cast<float>(psi_relaxed);
            }
        }
    }
    return std::sqrt(residual_l2 / static_cast<double>(m_W * m_H));
}

float DriftDiffusion::psiAt(int i, int j) const noexcept {
    if (i < 0 || j < 0 || i >= m_W || j >= m_H) return 0.0f;
    return m_psi[idx(i, j)];
}

float DriftDiffusion::bandShiftAt(int i, int j) const noexcept {
    // E_band(x,y) - E_band_ref = -q * psi(x,y).  q*psi gives J;
    // dividing by q gives eV; numerically -psi (V) is eV directly.
    return -psiAt(i, j);
}


// =============================================================================
// Electric field from psi gradient (central differences; one-sided at edges).
//   E_x = -(psi[i+1,j] - psi[i-1,j]) / (2 h)
//   E_y = -(psi[i,j+1] - psi[i,j-1]) / (2 h)
// Returned in V/cm.
// =============================================================================
float DriftDiffusion::electricFieldX(int i, int j) const noexcept {
    if (i < 0 || j < 0 || i >= m_W || j >= m_H) return 0.0f;
    const int ip = std::min(i + 1, m_W - 1);
    const int im = std::max(i - 1, 0);
    const float dx = static_cast<float>(ip - im) * m_cell_pitch_cm;
    if (dx == 0.0f) return 0.0f;
    return -(m_psi[idx(ip, j)] - m_psi[idx(im, j)]) / dx;
}

float DriftDiffusion::electricFieldY(int i, int j) const noexcept {
    if (i < 0 || j < 0 || i >= m_W || j >= m_H) return 0.0f;
    const int jp = std::min(j + 1, m_H - 1);
    const int jm = std::max(j - 1, 0);
    const float dy = static_cast<float>(jp - jm) * m_cell_pitch_cm;
    if (dy == 0.0f) return 0.0f;
    return -(m_psi[idx(i, jp)] - m_psi[idx(i, jm)]) / dy;
}

float DriftDiffusion::electricFieldMagAt(int i, int j) const noexcept {
    const float Ex = electricFieldX(i, j);
    const float Ey = electricFieldY(i, j);
    return std::sqrt(Ex * Ex + Ey * Ey);
}

float DriftDiffusion::peakElectricField() const noexcept {
    float peak = 0.0f;
    for (int j = 0; j < m_H; ++j) {
        for (int i = 0; i < m_W; ++i) {
            const float E = electricFieldMagAt(i, j);
            if (E > peak) peak = E;
        }
    }
    return peak;
}


// =============================================================================
// Spatial cuts [Phase 2: BandView]
// =============================================================================
bool DriftDiffusion::sampleHorizontalCut(
    int j_row,
    std::span<float> out_psi,
    std::span<float> out_Ec,
    std::span<float> out_Ev,
    double Ec0, double Ev0) const noexcept
{
    if (j_row < 0 || j_row >= m_H) return false;
    if (out_psi.size() < static_cast<std::size_t>(m_W)
        || out_Ec.size() < static_cast<std::size_t>(m_W)
        || out_Ev.size() < static_cast<std::size_t>(m_W)) return false;

    for (int i = 0; i < m_W; ++i) {
        const float psi = m_psi[idx(i, j_row)];
        out_psi[i] = psi;
        out_Ec [i] = static_cast<float>(Ec0) - psi;   // eV
        out_Ev [i] = static_cast<float>(Ev0) - psi;
    }
    return true;
}

bool DriftDiffusion::sampleEFieldCut(
    int j_row, std::span<float> out_E) const noexcept
{
    if (j_row < 0 || j_row >= m_H) return false;
    if (out_E.size() < static_cast<std::size_t>(m_W)) return false;
    for (int i = 0; i < m_W; ++i) {
        out_E[i] = electricFieldMagAt(i, j_row);
    }
    return true;
}

bool DriftDiffusion::sampleDopingCut(
    int j_row, std::span<float> out_net) const noexcept
{
    if (j_row < 0 || j_row >= m_H) return false;
    if (out_net.size() < static_cast<std::size_t>(m_W)) return false;
    for (int i = 0; i < m_W; ++i) {
        const auto k = idx(i, j_row);
        out_net[i] = m_Nd[k] - m_Na[k];
    }
    return true;
}


// =============================================================================
// Phase 3 -- Non-equilibrium / Gummel iteration
// =============================================================================
void DriftDiffusion::setAppliedBias(float V_a) noexcept {
    m_V_bias = std::clamp(V_a, -50.0f, 5.0f);
    // While AC is sweeping, the user dragging the bias slider should
    // re-centre the sine on the new DC point; otherwise the next step
    // would clobber the slider value with V_dc_base + sine.
    if (m_ac_enabled) m_V_dc_base = m_V_bias;
}

float DriftDiffusion::phiN(int i, int j) const noexcept {
    if (i < 0 || j < 0 || i >= m_W || j >= m_H) return 0.0f;
    return m_phi_n[idx(i, j)];
}
float DriftDiffusion::phiP(int i, int j) const noexcept {
    if (i < 0 || j < 0 || i >= m_W || j >= m_H) return 0.0f;
    return m_phi_p[idx(i, j)];
}

double DriftDiffusion::nDensityAt(int i, int j,
                                  double n_i, double V_T) const noexcept
{
    if (i < 0 || j < 0 || i >= m_W || j >= m_H) return 0.0;
    if (V_T <= 0.0) return 0.0;
    const auto k = idx(i, j);
    if (m_dens_seeded) return static_cast<double>(m_n_dens[k]);
    // Fallback (pre-seed): Boltzmann from phi.
    const double x = std::clamp(
        (static_cast<double>(m_psi[k]) - static_cast<double>(m_phi_n[k])) / V_T,
        -40.0, 40.0);
    return n_i * std::exp(x);
}
double DriftDiffusion::pDensityAt(int i, int j,
                                  double n_i, double V_T) const noexcept
{
    if (i < 0 || j < 0 || i >= m_W || j >= m_H) return 0.0;
    if (V_T <= 0.0) return 0.0;
    const auto k = idx(i, j);
    if (m_dens_seeded) return static_cast<double>(m_p_dens[k]);
    const double x = std::clamp(
        (static_cast<double>(m_phi_p[k]) - static_cast<double>(m_psi[k])) / V_T,
        -40.0, 40.0);
    return n_i * std::exp(x);
}


// =============================================================================
// Identify contacts on the leftmost / rightmost columns based on majority
// dopant.  Anode = column dominated by acceptors (Na > Nd); cathode =
// column dominated by donors. Single pass; cheap; runs once per Gummel
// invocation when m_contacts_dirty is set.
//
// Threshold: total dopant > 1e10 cm^-3 in the column (so an unintentional
// stray brush stamp does not turn a near-intrinsic edge into a contact).
// =============================================================================
namespace {

[[nodiscard]] bool columnIsContact(double sum) noexcept {
    return sum > 1.0e10;
}

} // namespace

void DriftDiffusion::applyContactBoundaries(double n_i, double V_T) noexcept
{
    if (n_i <= 0.0 || V_T <= 0.0) return;

    if (m_contacts_dirty) {
        std::fill(m_contact.begin(), m_contact.end(), 0u);

        // Aggregate Nd, Na for left and right outer columns.
        double sumLN = 0.0, sumLA = 0.0, sumRN = 0.0, sumRA = 0.0;
        for (int j = 0; j < m_H; ++j) {
            sumLN += m_Nd[idx(0,         j)];
            sumLA += m_Na[idx(0,         j)];
            sumRN += m_Nd[idx(m_W - 1,   j)];
            sumRA += m_Na[idx(m_W - 1,   j)];
        }

        const bool leftContact  = columnIsContact(sumLN + sumLA);
        const bool rightContact = columnIsContact(sumRN + sumRA);
        const bool leftAnode    = sumLA > sumLN;
        const bool rightAnode   = sumRA > sumRN;

        // Tag every cell of the qualifying outer column with the
        // appropriate role. m_contact: 1 = anode, 2 = cathode.
        if (leftContact) {
            const std::uint8_t tag = leftAnode ? 1u : 2u;
            for (int j = 0; j < m_H; ++j)
                m_contact[idx(0, j)] = tag;
        }
        if (rightContact) {
            const std::uint8_t tag = rightAnode ? 1u : 2u;
            for (int j = 0; j < m_H; ++j)
                m_contact[idx(m_W - 1, j)] = tag;
        }
        m_contacts_dirty = false;
    }

    // Apply Dirichlet values everywhere a contact tag is set. Anode is
    // at +V_a (forward = positive); cathode is grounded.
    //
    //   psi_bc   = V_metal +/- V_T ln(N_majority / n_i)   (Sze 2.7)
    //   n_bc     = N_d  (n-contact)  or  n_i^2 / N_a  (p-contact)
    //   p_bc     = n_i^2 / N_d (n-cont) or N_a (p-cont)
    //   phi_n_bc = phi_p_bc = V_metal     (no quasi-Fermi splitting at
    //                                       a metal interface)
    //
    // Phase 4: density Dirichlet is essential for the SG continuity
    // solver -- without it, n/p at the contacts decay to zero by
    // numerical dissipation and the device disconnects.
    const double n_i_sq = n_i * n_i;
    for (int j = 0; j < m_H; ++j) {
        for (int i = 0; i < m_W; ++i) {
            const std::size_t k = idx(i, j);
            if (m_contact[k] == 0u) continue;

            const double V_metal = (m_contact[k] == 1u)
                ? static_cast<double>(m_V_bias) : 0.0;

            const double Nd = std::max(double{m_Nd[k]}, 1.0);
            const double Na = std::max(double{m_Na[k]}, 1.0);

            const bool   isAnode = (m_contact[k] == 1u);
            const double psi_bc = isAnode
                ? V_metal - V_T * std::log(Na / n_i)                  // p-anode
                : V_metal + V_T * std::log(Nd / n_i);                 // n-cathode

            const double n_bc = isAnode ? n_i_sq / Na : Nd;
            const double p_bc = isAnode ? Na          : n_i_sq / Nd;

            m_psi   [k] = static_cast<float>(psi_bc);
            m_phi_n [k] = static_cast<float>(V_metal);
            m_phi_p [k] = static_cast<float>(V_metal);
            m_n_dens[k] = static_cast<float>(n_bc);
            m_p_dens[k] = static_cast<float>(p_bc);
        }
    }
}


// =============================================================================
// Density seeding & quasi-Fermi refresh (Phase 4)
// =============================================================================
void DriftDiffusion::seedDensitiesFromBoltzmann(
    double n_i, double V_T) noexcept
{
    if (n_i <= 0.0 || V_T <= 0.0) return;
    for (std::size_t k = 0; k < m_n_dens.size(); ++k) {
        const double xn = std::clamp(
            (static_cast<double>(m_psi[k]) - static_cast<double>(m_phi_n[k]))
            / V_T, -40.0, 40.0);
        const double xp = std::clamp(
            (static_cast<double>(m_phi_p[k]) - static_cast<double>(m_psi[k]))
            / V_T, -40.0, 40.0);
        m_n_dens[k] = static_cast<float>(n_i * std::exp(xn));
        m_p_dens[k] = static_cast<float>(n_i * std::exp(xp));
    }
    m_dens_seeded = true;
}

void DriftDiffusion::refreshQuasiFermiFromDensities(
    double n_i, double V_T) noexcept
{
    if (n_i <= 0.0 || V_T <= 0.0) return;
    constexpr float density_floor = 1.0e-30f;
    for (std::size_t k = 0; k < m_phi_n.size(); ++k) {
        const float n_C = std::max(m_n_dens[k], density_floor);
        const float p_C = std::max(m_p_dens[k], density_floor);
        m_phi_n[k] = static_cast<float>(
            static_cast<double>(m_psi[k])
            - V_T * std::log(static_cast<double>(n_C) / n_i));
        m_phi_p[k] = static_cast<float>(
            static_cast<double>(m_psi[k])
            + V_T * std::log(static_cast<double>(p_C) / n_i));
    }
}


// =============================================================================
// Scharfetter-Gummel continuity solvers              [Phase 4]
// -----------------------------------------------------------------------------
// SG flux between nodes (i,j) and (i+1,j) on uniform pitch h:
//
//   J_n(i+1/2) = (q V_T / h) mu_face [B(x) n_{i+1} - B(-x) n_i]
//   J_p(i+1/2) = -(q V_T / h) mu_face [B(-x) p_{i+1} - B(x) p_i]
//   x = (psi_{i+1} - psi_i) / V_T
//
// 5-point divergence + steady-state continuity gives the Gauss-Seidel
// update for n_C:
//
//   n_C^new = [alpha n_old_C + sum_k mu_k off_k n_k - h^2 (R - G) / V_T]
//           / [alpha + sum_k mu_k diag_k]
//
//   off_E   = B(x_E)        diag_E  = B(-x_E)            (electrons)
//   off_W   = B(-x_W)       diag_W  = B(x_W)
//   alpha   = h^2 / (V_T dt)   (zero -> steady-state limit)
//
// For holes, B(x) and B(-x) trade places (opposite drift sign).
//
// mu_face is evaluated *per face* via PhysicsEngine::localMobility*
// using the locally-averaged total ionised doping and the local field
// magnitude. This delivers velocity saturation inside depletion regions
// without any extra bookkeeping.
//
// The transient term (alpha != 0) is the Backward Euler form of
// (n_new - n_old)/dt = div(J_n)/q + (G - R), which is unconditionally
// stable -- dt may be picked freely from picoseconds to microseconds.
//
// Reference: Scharfetter-Gummel IEEE-ED 16 (1969) 64; Selberherr 6.2.
// =============================================================================
namespace {

// Local-mobility face evaluator. Pre-cached lambda factory keeps the
// hot-path branch-free; mat is captured by reference.
struct FaceMu {
    const material::Profile& mat;
    double T_lattice;
    double mu_bulk_n;
    double mu_bulk_p;

    // average face doping (cm^-3)
    [[nodiscard]] static double faceN(double Nd_a, double Nd_b,
                                      double Na_a, double Na_b) noexcept {
        return 0.5 * (Nd_a + Nd_b + Na_a + Na_b);
    }
};

} // namespace

void DriftDiffusion::solveContinuityElectron(
    double n_i, double V_T,
    double mu_n_bulk,
    const material::Profile& mat,
    double T_lattice,
    int sweeps, double omega) noexcept
{
    if (n_i <= 0.0 || V_T <= 0.0 || mu_n_bulk <= 0.0) return;
    if (!m_dens_seeded) seedDensitiesFromBoltzmann(n_i, V_T);

    const double h     = static_cast<double>(m_cell_pitch_cm);
    const double h2    = h * h;
    const double alpha = (m_transient_enabled && m_dt > 0.0)
                       ? h2 / (V_T * m_dt) : 0.0;

    for (int it = 0; it < sweeps; ++it) {
        for (int j = 0; j < m_H; ++j) {
            for (int i = 0; i < m_W; ++i) {
                const std::size_t k = idx(i, j);
                if (m_contact[k] != 0u) continue;

                const int ip = std::min(i + 1, m_W - 1);
                const int im = std::max(i - 1, 0);
                const int jp = std::min(j + 1, m_H - 1);
                const int jm = std::max(j - 1, 0);
                const std::size_t kE = idx(ip, j);
                const std::size_t kW = idx(im, j);
                const std::size_t kN = idx(i, jp);
                const std::size_t kS = idx(i, jm);

                const double psi_C = m_psi[k];
                const double xE = (m_psi[kE] - psi_C) / V_T;
                const double xW = (psi_C - m_psi[kW]) / V_T;
                const double xN = (m_psi[kN] - psi_C) / V_T;
                const double xS = (psi_C - m_psi[kS]) / V_T;

                // Per-face mobility (Matthiessen + Caughey-Thomas).
                // |E|_face from one-sided difference; field in V/cm.
                const double E_face_E = std::abs(m_psi[kE] - psi_C) / h;
                const double E_face_W = std::abs(psi_C - m_psi[kW]) / h;
                const double E_face_N = std::abs(m_psi[kN] - psi_C) / h;
                const double E_face_S = std::abs(psi_C - m_psi[kS]) / h;

                const double Nf_E = 0.5 * (m_Nd[k] + m_Nd[kE]
                                          + m_Na[k] + m_Na[kE]);
                const double Nf_W = 0.5 * (m_Nd[k] + m_Nd[kW]
                                          + m_Na[k] + m_Na[kW]);
                const double Nf_N = 0.5 * (m_Nd[k] + m_Nd[kN]
                                          + m_Na[k] + m_Na[kN]);
                const double Nf_S = 0.5 * (m_Nd[k] + m_Nd[kS]
                                          + m_Na[k] + m_Na[kS]);

                const double mu_E = (Nf_E > 1.0e10)
                    ? PhysicsEngine::localMobilityElectron(
                          mat, T_lattice, Nf_E, E_face_E)
                    : mu_n_bulk;
                const double mu_W = (Nf_W > 1.0e10)
                    ? PhysicsEngine::localMobilityElectron(
                          mat, T_lattice, Nf_W, E_face_W)
                    : mu_n_bulk;
                const double mu_N = (Nf_N > 1.0e10)
                    ? PhysicsEngine::localMobilityElectron(
                          mat, T_lattice, Nf_N, E_face_N)
                    : mu_n_bulk;
                const double mu_S = (Nf_S > 1.0e10)
                    ? PhysicsEngine::localMobilityElectron(
                          mat, T_lattice, Nf_S, E_face_S)
                    : mu_n_bulk;

                // SG Bernoulli coefficients.
                using PE = PhysicsEngine;
                const double bE  = PE::bernoulli( xE);
                const double bnE = PE::bernoulli(-xE);
                const double bW  = PE::bernoulli( xW);
                const double bnW = PE::bernoulli(-xW);
                const double bN  = PE::bernoulli( xN);
                const double bnN = PE::bernoulli(-xN);
                const double bS  = PE::bernoulli( xS);
                const double bnS = PE::bernoulli(-xS);

                // Electrons: off uses B(x_E)/B(-x_W)/B(x_N)/B(-x_S);
                //            diag uses B(-x_E)/B(x_W)/B(-x_N)/B(x_S).
                const double off_E = mu_E * bE;
                const double off_W = mu_W * bnW;
                const double off_N = mu_N * bN;
                const double off_S = mu_S * bnS;
                const double diag  = mu_E * bnE + mu_W * bW
                                   + mu_N * bnN + mu_S * bS;

                // Source: R_SRH + R_Aug - G_BTBT  [cm^-3 / s]
                const double n_C = m_n_dens[k];
                const double p_C = m_p_dens[k];
                const double R_srh = PE::recombSRH(
                    n_C, p_C, n_i, mat.tau_n, mat.tau_p);
                const double R_aug = PE::recombAuger(
                    n_C, p_C, n_i, mat.C_n_aug, mat.C_p_aug);
                const double Em    = electricFieldMagAt(i, j);
                const double G_bt  = PE::kaneBTBT(
                    Em, mat.A_kane, mat.B_kane, mat.btbt_isDirect);
                const double RmG   = R_srh + R_aug - G_bt;

                // n_C^new = [alpha n_old_C + sum off n_k - h^2 (R-G)/V_T]
                //        / [alpha + sum diag]
                const double rhs =
                    alpha * static_cast<double>(m_n_dens_old[k])
                  + off_E * m_n_dens[kE] + off_W * m_n_dens[kW]
                  + off_N * m_n_dens[kN] + off_S * m_n_dens[kS]
                  - h2 * RmG / V_T;
                const double den = alpha + diag;
                if (den <= 1.0e-30) continue;        // unreachable

                const double n_new = rhs / den;
                const double n_relaxed =
                    static_cast<double>(m_n_dens[k])
                  + omega * (n_new - static_cast<double>(m_n_dens[k]));

                // Floor at 0 -- B-coefficient products can occasionally
                // produce a negative excursion at the diffusion edge.
                m_n_dens[k] = static_cast<float>(std::max(n_relaxed, 0.0));
            }
        }
    }
}

void DriftDiffusion::solveContinuityHole(
    double n_i, double V_T,
    double mu_p_bulk,
    const material::Profile& mat,
    double T_lattice,
    int sweeps, double omega) noexcept
{
    if (n_i <= 0.0 || V_T <= 0.0 || mu_p_bulk <= 0.0) return;
    if (!m_dens_seeded) seedDensitiesFromBoltzmann(n_i, V_T);

    const double h     = static_cast<double>(m_cell_pitch_cm);
    const double h2    = h * h;
    const double alpha = (m_transient_enabled && m_dt > 0.0)
                       ? h2 / (V_T * m_dt) : 0.0;

    for (int it = 0; it < sweeps; ++it) {
        for (int j = 0; j < m_H; ++j) {
            for (int i = 0; i < m_W; ++i) {
                const std::size_t k = idx(i, j);
                if (m_contact[k] != 0u) continue;

                const int ip = std::min(i + 1, m_W - 1);
                const int im = std::max(i - 1, 0);
                const int jp = std::min(j + 1, m_H - 1);
                const int jm = std::max(j - 1, 0);
                const std::size_t kE = idx(ip, j);
                const std::size_t kW = idx(im, j);
                const std::size_t kN = idx(i, jp);
                const std::size_t kS = idx(i, jm);

                const double psi_C = m_psi[k];
                const double xE = (m_psi[kE] - psi_C) / V_T;
                const double xW = (psi_C - m_psi[kW]) / V_T;
                const double xN = (m_psi[kN] - psi_C) / V_T;
                const double xS = (psi_C - m_psi[kS]) / V_T;

                const double E_face_E = std::abs(m_psi[kE] - psi_C) / h;
                const double E_face_W = std::abs(psi_C - m_psi[kW]) / h;
                const double E_face_N = std::abs(m_psi[kN] - psi_C) / h;
                const double E_face_S = std::abs(psi_C - m_psi[kS]) / h;

                const double Nf_E = 0.5 * (m_Nd[k] + m_Nd[kE]
                                          + m_Na[k] + m_Na[kE]);
                const double Nf_W = 0.5 * (m_Nd[k] + m_Nd[kW]
                                          + m_Na[k] + m_Na[kW]);
                const double Nf_N = 0.5 * (m_Nd[k] + m_Nd[kN]
                                          + m_Na[k] + m_Na[kN]);
                const double Nf_S = 0.5 * (m_Nd[k] + m_Nd[kS]
                                          + m_Na[k] + m_Na[kS]);

                const double mu_E = (Nf_E > 1.0e10)
                    ? PhysicsEngine::localMobilityHole(
                          mat, T_lattice, Nf_E, E_face_E)
                    : mu_p_bulk;
                const double mu_W = (Nf_W > 1.0e10)
                    ? PhysicsEngine::localMobilityHole(
                          mat, T_lattice, Nf_W, E_face_W)
                    : mu_p_bulk;
                const double mu_N = (Nf_N > 1.0e10)
                    ? PhysicsEngine::localMobilityHole(
                          mat, T_lattice, Nf_N, E_face_N)
                    : mu_p_bulk;
                const double mu_S = (Nf_S > 1.0e10)
                    ? PhysicsEngine::localMobilityHole(
                          mat, T_lattice, Nf_S, E_face_S)
                    : mu_p_bulk;

                using PE = PhysicsEngine;
                const double bE  = PE::bernoulli( xE);
                const double bnE = PE::bernoulli(-xE);
                const double bW  = PE::bernoulli( xW);
                const double bnW = PE::bernoulli(-xW);
                const double bN  = PE::bernoulli( xN);
                const double bnN = PE::bernoulli(-xN);
                const double bS  = PE::bernoulli( xS);
                const double bnS = PE::bernoulli(-xS);

                // Holes: B(x) and B(-x) trade roles (opposite drift sign).
                const double off_E = mu_E * bnE;
                const double off_W = mu_W * bW;
                const double off_N = mu_N * bnN;
                const double off_S = mu_S * bS;
                const double diag  = mu_E * bE + mu_W * bnW
                                   + mu_N * bN + mu_S * bnS;

                const double n_C = m_n_dens[k];
                const double p_C = m_p_dens[k];
                const double R_srh = PE::recombSRH(
                    n_C, p_C, n_i, mat.tau_n, mat.tau_p);
                const double R_aug = PE::recombAuger(
                    n_C, p_C, n_i, mat.C_n_aug, mat.C_p_aug);
                const double Em    = electricFieldMagAt(i, j);
                const double G_bt  = PE::kaneBTBT(
                    Em, mat.A_kane, mat.B_kane, mat.btbt_isDirect);
                const double RmG   = R_srh + R_aug - G_bt;

                // p continuity: dp/dt = -div(J_p)/q + (G - R)
                // The SG sign on the source is the same as electrons
                // because the operator we derived already captures the
                // sign flip (B trade). Keep -h^2 (R-G)/V_T identical.
                const double rhs =
                    alpha * static_cast<double>(m_p_dens_old[k])
                  + off_E * m_p_dens[kE] + off_W * m_p_dens[kW]
                  + off_N * m_p_dens[kN] + off_S * m_p_dens[kS]
                  - h2 * RmG / V_T;
                const double den = alpha + diag;
                if (den <= 1.0e-30) continue;

                const double p_new = rhs / den;
                const double p_relaxed =
                    static_cast<double>(m_p_dens[k])
                  + omega * (p_new - static_cast<double>(m_p_dens[k]));
                m_p_dens[k] = static_cast<float>(std::max(p_relaxed, 0.0));
            }
        }
    }
}


// =============================================================================
// Gummel iteration
// -----------------------------------------------------------------------------
// Decoupled outer loop (Gummel, 1964): in each pass, we (1) lock the
// quasi-Fermi potentials and update psi from Poisson; (2) lock psi and
// update phi_n from electron continuity; (3) likewise phi_p from hole
// continuity. Each inner solver is Gauss-Seidel under-relaxation.
//
// Convergence indicator: max |delta psi| over the last outer iteration.
// (We accumulate the L2 residual in solvePoissonInner -- that already
// reports it, so we reuse it here as the return value.)
// =============================================================================
double DriftDiffusion::solveGummel(
    double n_i, double V_T, double epsilon_r,
    double mu_n_bulk, double mu_p_bulk,
    const material::Profile& mat,
    double T_lattice,
    int    outer_iters,
    int    poisson_inner,
    int    continuity_inner,
    double omega_psi,
    double omega_phi) noexcept
{
    if (outer_iters <= 0) return 0.0;

    // Bootstrap densities on first call (or after clearDoping).
    if (!m_dens_seeded) {
        applyContactBoundaries(n_i, V_T);
        seedDensitiesFromBoltzmann(n_i, V_T);
    }

    double residual = 0.0;
    for (int outer = 0; outer < outer_iters; ++outer) {
        applyContactBoundaries(n_i, V_T);
        residual = solvePoissonInner(n_i, V_T, epsilon_r,
                                     poisson_inner, omega_psi);
        applyContactBoundaries(n_i, V_T);    // refresh ψ after sweep

        solveContinuityElectron(n_i, V_T, mu_n_bulk, mat, T_lattice,
                                continuity_inner, omega_phi);
        solveContinuityHole    (n_i, V_T, mu_p_bulk, mat, T_lattice,
                                continuity_inner, omega_phi);

        // Sync visualisation / Poisson view: phi_n,p must reflect the
        // SG-iterated densities so the next Poisson sweep sees the same
        // n,p values. Poisson reads phi -> exp -> n; with this sync the
        // Boltzmann formula and the SG densities agree at this iterate.
        refreshQuasiFermiFromDensities(n_i, V_T);
    }
    return residual;
}


// =============================================================================
// Quasi-Fermi cut sampler (BandView)
// =============================================================================
bool DriftDiffusion::sampleQuasiFermiCut(
    int j_row,
    std::span<float> out_phi_n,
    std::span<float> out_phi_p) const noexcept
{
    if (j_row < 0 || j_row >= m_H) return false;
    if (out_phi_n.size() < static_cast<std::size_t>(m_W)
        || out_phi_p.size() < static_cast<std::size_t>(m_W)) return false;
    for (int i = 0; i < m_W; ++i) {
        const auto k = idx(i, j_row);
        out_phi_n[i] = m_phi_n[k];
        out_phi_p[i] = m_phi_p[k];
    }
    return true;
}


// =============================================================================
// Terminal current density (mid-column proxy)
// -----------------------------------------------------------------------------
// J_n = -q mu_n n d phi_n / dx
// J_p = -q mu_p p d phi_p / dx     (same sign convention)
// Total J = J_n + J_p, averaged along i = W/2 column.
// Returns A/cm^2.
// =============================================================================
double DriftDiffusion::terminalCurrentDensity(
    double n_i, double V_T,
    double mu_n, double mu_p) const noexcept
{
    if (n_i <= 0.0 || V_T <= 0.0 || m_W < 3 || m_H < 1) return 0.0;
    constexpr double q = 1.602176634e-19;
    const double h    = static_cast<double>(m_cell_pitch_cm);

    // Accumulate the SG flux directly on a mid-column face -- this is
    // the quantity the SG continuity iterates already conserve, so the
    // result is consistent with the iterated state. (Previous central-
    // difference of phi worked but mismatched the conservation form.)
    const int i_mid = m_W / 2;
    if (i_mid < 1 || i_mid >= m_W) return 0.0;

    double J_sum = 0.0;
    for (int j = 0; j < m_H; ++j) {
        const auto kL = idx(i_mid,     j);
        const auto kR = idx(i_mid + 1, j);

        const double x = (m_psi[kR] - m_psi[kL]) / V_T;
        const double bp = PhysicsEngine::bernoulli( x);
        const double bn = PhysicsEngine::bernoulli(-x);

        const double J_n = (q * mu_n * V_T / h)
            * (bp * m_n_dens[kR] - bn * m_n_dens[kL]);
        const double J_p = -(q * mu_p * V_T / h)
            * (bn * m_p_dens[kR] - bp * m_p_dens[kL]);

        J_sum += J_n + J_p;
    }
    return J_sum / static_cast<double>(m_H);
}


// =============================================================================
// Transient (Backward Euler) and AC probe              [Phase 4]
// =============================================================================
void DriftDiffusion::setTimeStep(double dt_s) noexcept {
    // Clamp to a reasonable range so the BE transient term doesn't get
    // pathological. Picoseconds at the low end (tunnel-current scales);
    // milliseconds at the high end (long thermal-diffusion-style runs).
    m_dt = std::clamp(dt_s, 1.0e-15, 1.0e-3);
}

void DriftDiffusion::setACProbe(bool enabled,
                                double freq_Hz, double amp_V) noexcept
{
    m_ac_enabled = enabled;
    m_ac_freq    = std::clamp(freq_Hz, 1.0, 1.0e12);
    m_ac_amp     = std::clamp(amp_V,   0.0, 1.0);
    if (enabled) m_V_dc_base = m_V_bias;
}

double DriftDiffusion::stepTransient(
    double n_i, double V_T, double epsilon_r,
    const material::Profile& mat,
    double T_lattice,
    int    outer_iters,
    int    poisson_inner,
    int    continuity_inner,
    double omega_psi,
    double omega_phi) noexcept
{
    if (!m_transient_enabled || m_dt <= 0.0) return 0.0;

    // Seed densities first if needed; otherwise copying zero -> n_old
    // would feed a zero into the BE diagonal and crush the first step.
    if (!m_dens_seeded) {
        applyContactBoundaries(n_i, V_T);
        seedDensitiesFromBoltzmann(n_i, V_T);
    }

    // 1) Save the current state into the "old" buffers for the BE term.
    std::copy(m_n_dens.begin(), m_n_dens.end(), m_n_dens_old.begin());
    std::copy(m_p_dens.begin(), m_p_dens.end(), m_p_dens_old.begin());

    // 2) Advance simulation time and apply the AC probe to V_a.
    m_sim_time += m_dt;
    if (m_ac_enabled) {
        constexpr double TWO_PI = 6.283185307179586;
        const double sineV = m_ac_amp
            * std::sin(TWO_PI * m_ac_freq * m_sim_time);
        m_V_bias = static_cast<float>(m_V_dc_base + sineV);
    }

    // 3) Bulk mobilities for the Matthiessen fallback at intrinsic faces.
    using PE = PhysicsEngine;
    const double mu_n_bulk = PE::matthiessenMobilityElectron(mat, T_lattice, 0.0);
    const double mu_p_bulk = PE::matthiessenMobilityHole    (mat, T_lattice, 0.0);

    // 4) Run one Gummel iteration with the BE transient term active.
    return solveGummel(n_i, V_T, epsilon_r,
                       mu_n_bulk, mu_p_bulk, mat, T_lattice,
                       outer_iters, poisson_inner, continuity_inner,
                       omega_psi, omega_phi);
}


// =============================================================================
// Small-signal conductance via DC perturbation          [Phase 4 bonus]
// -----------------------------------------------------------------------------
// Symmetric 2-point finite difference around the current bias:
//
//     G ~ ( J(V_a + dV) - J(V_a - dV) ) / (2 dV)        [S/cm^2]
//
// Cheap and robust. For C, the engine offers PhysicsEngine::depletion-
// CapacitanceFlat / diffusionCapacitance using the painted Nd, Na and
// the current forward I -- those run from the UI directly because they
// don't need a perturbation sweep.
// =============================================================================
double DriftDiffusion::smallSignalConductance(
    double n_i, double V_T, double epsilon_r,
    const material::Profile& mat,
    double T_lattice,
    double dV) noexcept
{
    if (dV <= 0.0) return 0.0;
    using PE = PhysicsEngine;
    const double mu_n = PE::matthiessenMobilityElectron(mat, T_lattice, 0.0);
    const double mu_p = PE::matthiessenMobilityHole    (mat, T_lattice, 0.0);

    // Save current state -- the perturbation runs trash the densities,
    // so we restore everything at the end.
    const float V_save = m_V_bias;
    static thread_local std::vector<float> save_n, save_p, save_psi,
                                            save_phi_n, save_phi_p;
    save_n    .resize(m_n_dens.size());   std::copy(m_n_dens.begin(),  m_n_dens.end(),  save_n.begin());
    save_p    .resize(m_p_dens.size());   std::copy(m_p_dens.begin(),  m_p_dens.end(),  save_p.begin());
    save_psi  .resize(m_psi.size());      std::copy(m_psi.begin(),     m_psi.end(),     save_psi.begin());
    save_phi_n.resize(m_phi_n.size());    std::copy(m_phi_n.begin(),   m_phi_n.end(),   save_phi_n.begin());
    save_phi_p.resize(m_phi_p.size());    std::copy(m_phi_p.begin(),   m_phi_p.end(),   save_phi_p.begin());

    constexpr int OUTER  = 12;
    constexpr int POISS  = 30;
    constexpr int CONT   = 18;

    // J at V + dV
    setAppliedBias(V_save + static_cast<float>(dV));
    solveGummel(n_i, V_T, epsilon_r, mu_n, mu_p, mat, T_lattice,
                OUTER, POISS, CONT, 0.85, 1.0);
    const double J_plus  = terminalCurrentDensity(n_i, V_T, mu_n, mu_p);

    // restore mid-point
    std::copy(save_n.begin(),     save_n.end(),     m_n_dens.begin());
    std::copy(save_p.begin(),     save_p.end(),     m_p_dens.begin());
    std::copy(save_psi.begin(),   save_psi.end(),   m_psi.begin());
    std::copy(save_phi_n.begin(), save_phi_n.end(), m_phi_n.begin());
    std::copy(save_phi_p.begin(), save_phi_p.end(), m_phi_p.begin());

    // J at V - dV
    setAppliedBias(V_save - static_cast<float>(dV));
    solveGummel(n_i, V_T, epsilon_r, mu_n, mu_p, mat, T_lattice,
                OUTER, POISS, CONT, 0.85, 1.0);
    const double J_minus = terminalCurrentDensity(n_i, V_T, mu_n, mu_p);

    // Restore.
    std::copy(save_n.begin(),     save_n.end(),     m_n_dens.begin());
    std::copy(save_p.begin(),     save_p.end(),     m_p_dens.begin());
    std::copy(save_psi.begin(),   save_psi.end(),   m_psi.begin());
    std::copy(save_phi_n.begin(), save_phi_n.end(), m_phi_n.begin());
    std::copy(save_phi_p.begin(), save_phi_p.end(), m_phi_p.begin());
    setAppliedBias(V_save);

    return (J_plus - J_minus) / (2.0 * dV);
}
