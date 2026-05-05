// =============================================================================
// DriftDiffusion.cpp
//
//   Author : dex / cemfly-april2026
//   License: MIT
// =============================================================================

#include "DriftDiffusion.hpp"

#include <algorithm>
#include <cmath>


// =============================================================================
// Construction -- pre-allocate every grid buffer once and never again.
// =============================================================================
DriftDiffusion::DriftDiffusion(int gridW, int gridH)
    : m_W(gridW), m_H(gridH)
    , m_n     (static_cast<std::size_t>(gridW * gridH), 0.0f)
    , m_n_next(static_cast<std::size_t>(gridW * gridH), 0.0f)
    , m_G     (static_cast<std::size_t>(gridW * gridH), 0.0f)
    , m_T     (static_cast<std::size_t>(gridW * gridH), 300.0f)
    , m_T_next(static_cast<std::size_t>(gridW * gridH), 300.0f)
    , m_region(static_cast<std::size_t>(gridW * gridH),
               static_cast<std::uint8_t>(CellRegion::Bulk))
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
// Source helpers (unchanged from Phase 5)
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
    rebuildRegionMap();
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
    constexpr float k_B_T_300 = 0.02585f;     // kT at 300 K [eV]
    const float thermalArg = std::clamp(m_VBE / k_B_T_300, 0.0f, 25.0f);
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
