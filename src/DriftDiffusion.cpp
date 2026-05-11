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
    // Phase 5: per-cell material id (defaults to "reference" = id 0xFF)
    , m_mat         (static_cast<std::size_t>(gridW * gridH), 0xFFu)
    , m_ni_local    (static_cast<std::size_t>(gridW * gridH), 1.0e10f)  // Si@300K placeholder
    , m_eps_r_local (static_cast<std::size_t>(gridW * gridH), 11.7f)
    , m_delta_n     (static_cast<std::size_t>(gridW * gridH), 0.0f)
    , m_delta_p     (static_cast<std::size_t>(gridW * gridH), 0.0f)
    // Phase 5: Wachutka thermal extras
    , m_T_old      (static_cast<std::size_t>(gridW * gridH), 300.0f)
    , m_H_gen      (static_cast<std::size_t>(gridW * gridH), 0.0f)
    // Phase 7: Newton-Krylov scratch (3 unknowns per cell -- doubles for
    // the linear-solver inner products to stay numerically stable).
    , m_nk_x      (static_cast<std::size_t>(3 * gridW * gridH), 0.0)
    , m_nk_dx     (static_cast<std::size_t>(3 * gridW * gridH), 0.0)
    , m_nk_F      (static_cast<std::size_t>(3 * gridW * gridH), 0.0)
    , m_nk_F_pert (static_cast<std::size_t>(3 * gridW * gridH), 0.0)
    , m_nk_x_pert (static_cast<std::size_t>(3 * gridW * gridH), 0.0)
    , m_bicg_r    (static_cast<std::size_t>(3 * gridW * gridH), 0.0)
    , m_bicg_rh   (static_cast<std::size_t>(3 * gridW * gridH), 0.0)
    , m_bicg_p    (static_cast<std::size_t>(3 * gridW * gridH), 0.0)
    , m_bicg_v    (static_cast<std::size_t>(3 * gridW * gridH), 0.0)
    , m_bicg_s    (static_cast<std::size_t>(3 * gridW * gridH), 0.0)
    , m_bicg_t    (static_cast<std::size_t>(3 * gridW * gridH), 0.0)
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

    // Phase 5: Identify which material::Kind this profile corresponds
    // to so painted reference cells (id 0xFF) resolve to the right
    // ChiSet during the heterojunction Poisson / Wachutka stages.
    if      (&mat == &material::byKind(material::Kind::Silicon))
        m_ref_mat = material::Kind::Silicon;
    else if (&mat == &material::byKind(material::Kind::GaAs))
        m_ref_mat = material::Kind::GaAs;
    else if (&mat == &material::byKind(material::Kind::Germanium))
        m_ref_mat = material::Kind::Germanium;
    // else: leave the previous reference -- caller wired a custom Profile,
    // and changing it would silently break the heterojunction map.
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
    std::fill(m_T_old.begin(),  m_T_old.end(),  m_T_ambient);
    std::fill(m_H_gen.begin(),  m_H_gen.end(),  0.0f);
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

void DriftDiffusion::setDopingAt(int i, int j, double Nd, double Na) noexcept {
    if (i < 0 || j < 0 || i >= m_W || j >= m_H) return;
    const std::size_t k = idx(i, j);
    m_Nd[k] = static_cast<float>(std::max(0.0, Nd));
    m_Na[k] = static_cast<float>(std::max(0.0, Na));
    // Region tag follows majority dopant; both-zero -> bulk.
    if      (Nd > Na && Nd > 0.0)
        m_region[k] = static_cast<std::uint8_t>(CellRegion::NDoped);
    else if (Na > Nd && Na > 0.0)
        m_region[k] = static_cast<std::uint8_t>(CellRegion::PDoped);
    else
        m_region[k] = static_cast<std::uint8_t>(CellRegion::Bulk);
    m_contacts_dirty = true;
}


// =============================================================================
// Heterojunction painter [Phase 5]
// -----------------------------------------------------------------------------
// Per-cell material id storage:
//
//   0xFF = "use the engine's reference material" (default everywhere)
//   0   = Silicon
//   1   = GaAs
//   2   = Germanium
//
// Storing 0xFF for unpainted cells lets us repaint the canvas after the
// engine swaps materials without touching every cell -- the resolveAt()
// helper looks up the reference profile only when the id is 0xFF.
// =============================================================================
void DriftDiffusion::setReferenceMaterial(material::Kind k) noexcept {
    m_ref_mat = k;
    // Painted cells keep their explicit id; only unpainted (0xFF) cells
    // follow the new reference -- which is what the user expects when
    // they swap the global Material dropdown after painting an HBT.
}

void DriftDiffusion::paintMaterialBrush(float u, float v,
                                        material::Kind k,
                                        int radius_cells) noexcept
{
    u = std::clamp(u, 0.0f, 1.0f);
    v = std::clamp(v, 0.0f, 1.0f);
    const int ic = static_cast<int>(u * (m_W - 1));
    const int jc = static_cast<int>(v * (m_H - 1));
    const int r2 = radius_cells * radius_cells;

    const int i0 = std::max(0,        ic - radius_cells);
    const int i1 = std::min(m_W - 1,  ic + radius_cells);
    const int j0 = std::max(0,        jc - radius_cells);
    const int j1 = std::min(m_H - 1,  jc + radius_cells);

    const auto id = static_cast<std::uint8_t>(k);
    for (int j = j0; j <= j1; ++j) {
        for (int i = i0; i <= i1; ++i) {
            const int di = i - ic;
            const int dj = j - jc;
            if (di * di + dj * dj > r2) continue;
            m_mat[idx(i, j)] = id;
        }
    }
    m_contacts_dirty = true;   // contact n_i depends on local material
}

void DriftDiffusion::clearMaterials() noexcept {
    std::fill(m_mat.begin(), m_mat.end(), 0xFFu);
    m_contacts_dirty = true;
}

material::Kind DriftDiffusion::materialAt(int i, int j) const noexcept {
    if (i < 0 || j < 0 || i >= m_W || j >= m_H) return m_ref_mat;
    const auto id = m_mat[idx(i, j)];
    if (id == 0xFFu) return m_ref_mat;
    if (id < static_cast<std::uint8_t>(material::kCount))
        return static_cast<material::Kind>(id);
    return m_ref_mat;
}

const material::Profile& DriftDiffusion::profileAt(int i, int j) const noexcept {
    return material::byKind(materialAt(i, j));
}

double DriftDiffusion::localBandgapAt(int i, int j) const noexcept {
    const auto& mat = profileAt(i, j);
    const float T   = temperatureFieldAt(i, j);
    return PhysicsEngine::bandgapAt(mat, static_cast<double>(T));
}

double DriftDiffusion::localIntrinsicAt(int i, int j) const noexcept {
    const auto& mat = profileAt(i, j);
    const float T   = temperatureFieldAt(i, j);
    return PhysicsEngine::intrinsicCarrierAt(mat, static_cast<double>(T));
}

double DriftDiffusion::localElectronAffinityAt(int i, int j) const noexcept {
    return profileAt(i, j).chi;
}

float DriftDiffusion::temperatureFieldAt(int i, int j) const noexcept {
    if (i < 0 || j < 0 || i >= m_W || j >= m_H) return m_T_ambient;
    return m_T[idx(i, j)];
}


// =============================================================================
// Heterojunction parameter cache refresh
// -----------------------------------------------------------------------------
// Pulls per-cell material parameters once, so the Poisson and continuity
// inner sweeps only ever touch float buffers (no virtual lookups, no
// std::pow / std::log10 in the hot path).
//
// Cost: O(W*H) per call; called at the start of every public solver entry
// point (solvePoisson, solveGummel, stepTransient, solveHeatEquation).
// =============================================================================
void DriftDiffusion::refreshLocalParamCache(double V_T) noexcept {
    (void)V_T;  // currently unused; kept for API stability vs future BGN

    const auto& mat_ref = material::byKind(m_ref_mat);
    const double chi_ref = mat_ref.chi;
    // Reference Eg evaluated at the device's *ambient* temperature. We do
    // not let Eg_ref drift with local hotspots -- the Anderson chi-shift
    // is a metallurgical property of the interface, not a hot-spot effect.
    const double Eg_ref  = PhysicsEngine::bandgapAt(
        mat_ref, static_cast<double>(m_T_ambient));

    for (int j = 0; j < m_H; ++j) {
        for (int i = 0; i < m_W; ++i) {
            const std::size_t k = idx(i, j);
            const auto& mat = profileAt(i, j);
            const double T  = static_cast<double>(m_T[k]);

            // Local n_i with hot-spot feedback: lattice T enters via both
            // bandgap narrowing (Varshni) and the T^(3/2) DOS scaling.
            // This is what closes the thermal-runaway loop.
            m_ni_local[k]    = static_cast<float>(
                PhysicsEngine::intrinsicCarrierAt(mat, T));

            // Local epsilon_r for the Poisson harmonic-mean stencil.
            m_eps_r_local[k] = static_cast<float>(mat.epsilon_r);

            // Anderson chi-shifts.  Eg used for delta_p uses ambient T
            // (same reasoning as Eg_ref above) so the sole T dependence
            // sits in n_i and the continuity-side mobility evaluator.
            const double Eg_C  = PhysicsEngine::bandgapAt(
                mat, static_cast<double>(m_T_ambient));
            const double dn = mat.chi - chi_ref;
            const double dp = (mat.chi - chi_ref) + (Eg_C - Eg_ref);
            m_delta_n[k] = static_cast<float>(dn);
            m_delta_p[k] = static_cast<float>(dp);
        }
    }
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
    // Refresh the heterojunction cache & contact BCs first so stale
    // state can't leak in.
    refreshLocalParamCache(V_T);
    applyContactBoundaries(n_i_cm3, V_T);
    return solvePoissonInner(n_i_cm3, V_T, epsilon_r, iterations, omega);
}

// =============================================================================
// Heterojunction-aware Poisson inner sweep                       [Phase 5]
// -----------------------------------------------------------------------------
// Finite-volume 5-point stencil with variable epsilon (harmonic mean):
//
//   eps_e (psi_E - psi_C) + eps_w (psi_W - psi_C)
// + eps_n (psi_N - psi_C) + eps_s (psi_S - psi_C)
//   = -h^2 q (p - n + N_d - N_a)
//
// solving for psi_C (Gauss-Seidel under-relaxation):
//
//   psi_C^new = [ eps_e psi_E + eps_w psi_W + eps_n psi_N + eps_s psi_S
//               + (q h^2 / eps_0) rho ] / (eps_e + eps_w + eps_n + eps_s)
//
// where the *relative* permittivities sit inside the harmonic mean and
// q*h^2/eps_0 carries the absolute scaling (volts).  In a homogeneous
// reference material, eps_e = eps_w = eps_n = eps_s = eps_r and the
// formula collapses to the Phase-1 expression with rho_coef = q h^2 / eps_s.
//
// The carrier densities use the Anderson chi-shifts m_delta_n / m_delta_p:
//
//   n(x) = n_i_ref * exp((psi + delta_n - phi_n) / V_T)
//   p(x) = n_i_ref * exp((phi_p - psi - delta_p) / V_T)
//
// so that np = n_i_ref^2 exp((delta_n - delta_p)/V_T) = n_i(x)^2 in
// equilibrium (mass action holds locally even though the prefactor is
// uniform).  See refreshLocalParamCache() for the definitions of
// delta_n / delta_p; in the homogeneous limit both are identically zero.
// References: Selberherr Sec. 5; Pierret App. B; Anderson 1962.
// =============================================================================
double DriftDiffusion::solvePoissonInner(double n_i_cm3,
                                         double V_T,
                                         double epsilon_r,
                                         int    sweeps,
                                         double omega) noexcept
{
    (void)epsilon_r;  // now per-cell; argument retained for API stability
    constexpr double eps_0_local = 8.854187817e-14;     // F/cm
    constexpr double q           = 1.602176634e-19;
    const double h               = static_cast<double>(m_cell_pitch_cm);
    const double rho_coef_pre    = q * h * h / eps_0_local;   // [V * cm^3]

    if (V_T <= 0.0) return 0.0;

    double residual_l2 = 0.0;

    for (int it = 0; it < sweeps; ++it) {
        residual_l2 = 0.0;

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

                // Per-cell relative permittivities and harmonic-mean
                // face values (Phase 5: heterojunction support).
                const double er_C = m_eps_r_local[k];
                const double er_E = m_eps_r_local[kE];
                const double er_W = m_eps_r_local[kW];
                const double er_N = m_eps_r_local[kN];
                const double er_S = m_eps_r_local[kS];

                const double ef_E = 2.0 * er_C * er_E / (er_C + er_E);
                const double ef_W = 2.0 * er_C * er_W / (er_C + er_W);
                const double ef_N = 2.0 * er_C * er_N / (er_C + er_N);
                const double ef_S = 2.0 * er_C * er_S / (er_C + er_S);

                // Anderson chi-shifted Boltzmann densities at the cell.
                // (n_i_cm3 used here as the reference n_i; see header.)
                const double dn   = m_delta_n[k];
                const double dp   = m_delta_p[k];
                const double psi_C  = m_psi[k];
                const double phin_C = m_phi_n[k];
                const double phip_C = m_phi_p[k];

                const double xn = std::clamp(
                    (psi_C + dn - phin_C) / V_T, -40.0, 40.0);
                const double xp = std::clamp(
                    (phip_C - psi_C - dp) / V_T, -40.0, 40.0);
                const double n_C = n_i_cm3 * std::exp(xn);
                const double p_C = n_i_cm3 * std::exp(xp);

                const double rho =
                    p_C - n_C
                    + static_cast<double>(m_Nd[k])
                    - static_cast<double>(m_Na[k]);

                const double psi_E = m_psi[kE];
                const double psi_W = m_psi[kW];
                const double psi_N = m_psi[kN];
                const double psi_S = m_psi[kS];

                const double eps_sum = ef_E + ef_W + ef_N + ef_S;
                const double psi_new = (
                      ef_E * psi_E + ef_W * psi_W
                    + ef_N * psi_N + ef_S * psi_S
                    + rho_coef_pre * rho) / eps_sum;

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

    // Phase 5: heterojunction-aware band edges.
    //
    //   E_c(x) = Ec0 - psi(x) - (chi(x) - chi_ref)
    //   E_v(x) = E_c(x) - E_g(x, T_local)
    //
    // The (chi - chi_ref) term is the Anderson conduction-band offset;
    // the local E_g(T) introduces both the gap-narrowing hot-spot effect
    // and the heterojunction valence-band step (since dEv = dEg - dEc).
    // In the homogeneous-reference limit chi == chi_ref and E_g == Ec0,
    // and the formula collapses to the Phase-2 expression bit-for-bit.
    (void)Ev0;
    const auto& mat_ref = material::byKind(m_ref_mat);
    const double chi_ref = mat_ref.chi;

    for (int i = 0; i < m_W; ++i) {
        const auto kij = idx(i, j_row);
        const float psi = m_psi[kij];
        out_psi[i] = psi;

        const auto& mat = profileAt(i, j_row);
        const double T_C  = static_cast<double>(m_T[kij]);
        const double Eg_C = PhysicsEngine::bandgapAt(mat, T_C);
        const double dn   = mat.chi - chi_ref;          // chi shift [eV]

        const float Ec = static_cast<float>(Ec0 - psi - dn);
        out_Ec[i] = Ec;
        out_Ev[i] = static_cast<float>(Ec - Eg_C);
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
    // m_dt may have been updated since the previous call (transient
    // sweeps). Read once per sweep set.
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

                // Phase 5: SG argument uses the *electron drift potential*
                // eta_n = psi + delta_n where delta_n = (chi - chi_ref)/q.
                // In the homogeneous-reference limit delta_n = 0 and we
                // recover the Phase-4 formulation bit-for-bit.
                const double dn_C = m_delta_n[k];
                const double psi_eta_C = m_psi[k]  + dn_C;
                const double psi_eta_E = m_psi[kE] + m_delta_n[kE];
                const double psi_eta_W = m_psi[kW] + m_delta_n[kW];
                const double psi_eta_N = m_psi[kN] + m_delta_n[kN];
                const double psi_eta_S = m_psi[kS] + m_delta_n[kS];

                const double xE = (psi_eta_E - psi_eta_C) / V_T;
                const double xW = (psi_eta_C - psi_eta_W) / V_T;
                const double xN = (psi_eta_N - psi_eta_C) / V_T;
                const double xS = (psi_eta_C - psi_eta_S) / V_T;

                // Per-face mobility (Matthiessen + Caughey-Thomas) using
                // the *local* material profile at each face.
                const double E_face_E = std::abs(psi_eta_E - psi_eta_C) / h;
                const double E_face_W = std::abs(psi_eta_C - psi_eta_W) / h;
                const double E_face_N = std::abs(psi_eta_N - psi_eta_C) / h;
                const double E_face_S = std::abs(psi_eta_C - psi_eta_S) / h;

                const double Nf_E = 0.5 * (m_Nd[k] + m_Nd[kE]
                                          + m_Na[k] + m_Na[kE]);
                const double Nf_W = 0.5 * (m_Nd[k] + m_Nd[kW]
                                          + m_Na[k] + m_Na[kW]);
                const double Nf_N = 0.5 * (m_Nd[k] + m_Nd[kN]
                                          + m_Na[k] + m_Na[kN]);
                const double Nf_S = 0.5 * (m_Nd[k] + m_Nd[kS]
                                          + m_Na[k] + m_Na[kS]);

                // Local lattice T per cell (Wachutka feedback) and per-
                // cell material for mobility evaluation.
                const auto& mat_C  = profileAt(i,  j);
                const double T_C   = static_cast<double>(m_T[k]);
                const auto& mat_E_ = profileAt(ip, j);
                const auto& mat_W_ = profileAt(im, j);
                const auto& mat_N_ = profileAt(i,  jp);
                const auto& mat_S_ = profileAt(i,  jm);

                const double mu_E = (Nf_E > 1.0e10)
                    ? PhysicsEngine::localMobilityElectron(
                          mat_E_, T_C, Nf_E, E_face_E)
                    : mu_n_bulk;
                const double mu_W = (Nf_W > 1.0e10)
                    ? PhysicsEngine::localMobilityElectron(
                          mat_W_, T_C, Nf_W, E_face_W)
                    : mu_n_bulk;
                const double mu_N = (Nf_N > 1.0e10)
                    ? PhysicsEngine::localMobilityElectron(
                          mat_N_, T_C, Nf_N, E_face_N)
                    : mu_n_bulk;
                const double mu_S = (Nf_S > 1.0e10)
                    ? PhysicsEngine::localMobilityElectron(
                          mat_S_, T_C, Nf_S, E_face_S)
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

                // Source:  U_net (= R_SRH + R_Aug + R_rad)  -  G_BTBT.
                // Local n_i depends on local Eg + lattice T (already
                // baked into m_ni_local). Radiative coefficient B_rad
                // varies by material (5 orders Si vs GaAs) -- read off
                // mat_C so per-cell heterostructures stay accurate.
                const double n_C = m_n_dens[k];
                const double p_C = m_p_dens[k];
                const double n_i_C = static_cast<double>(m_ni_local[k]);
                const double U     = PE::netRecombination(
                    n_C, p_C, n_i_C, mat_C);
                const double Em    = electricFieldMagAt(i, j);
                const double G_bt  = PE::kaneBTBT(
                    Em, mat_C.A_kane, mat_C.B_kane, mat_C.btbt_isDirect);
                const double RmG   = U - G_bt;
                (void)mat; (void)T_lattice;

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

                // Phase 5: SG argument uses the *hole drift potential*
                // eta_p = psi + delta_p where delta_p includes both the
                // chi step and the Eg step (-> -E_v).  This is what
                // makes a wide-gap emitter reject back-injected holes.
                const double psi_eta_C = m_psi[k]  + m_delta_p[k];
                const double psi_eta_E = m_psi[kE] + m_delta_p[kE];
                const double psi_eta_W = m_psi[kW] + m_delta_p[kW];
                const double psi_eta_N = m_psi[kN] + m_delta_p[kN];
                const double psi_eta_S = m_psi[kS] + m_delta_p[kS];

                const double xE = (psi_eta_E - psi_eta_C) / V_T;
                const double xW = (psi_eta_C - psi_eta_W) / V_T;
                const double xN = (psi_eta_N - psi_eta_C) / V_T;
                const double xS = (psi_eta_C - psi_eta_S) / V_T;

                const double E_face_E = std::abs(psi_eta_E - psi_eta_C) / h;
                const double E_face_W = std::abs(psi_eta_C - psi_eta_W) / h;
                const double E_face_N = std::abs(psi_eta_N - psi_eta_C) / h;
                const double E_face_S = std::abs(psi_eta_C - psi_eta_S) / h;

                const double Nf_E = 0.5 * (m_Nd[k] + m_Nd[kE]
                                          + m_Na[k] + m_Na[kE]);
                const double Nf_W = 0.5 * (m_Nd[k] + m_Nd[kW]
                                          + m_Na[k] + m_Na[kW]);
                const double Nf_N = 0.5 * (m_Nd[k] + m_Nd[kN]
                                          + m_Na[k] + m_Na[kN]);
                const double Nf_S = 0.5 * (m_Nd[k] + m_Nd[kS]
                                          + m_Na[k] + m_Na[kS]);

                const auto& mat_C  = profileAt(i,  j);
                const double T_C   = static_cast<double>(m_T[k]);
                const auto& mat_E_ = profileAt(ip, j);
                const auto& mat_W_ = profileAt(im, j);
                const auto& mat_N_ = profileAt(i,  jp);
                const auto& mat_S_ = profileAt(i,  jm);

                const double mu_E = (Nf_E > 1.0e10)
                    ? PhysicsEngine::localMobilityHole(
                          mat_E_, T_C, Nf_E, E_face_E)
                    : mu_p_bulk;
                const double mu_W = (Nf_W > 1.0e10)
                    ? PhysicsEngine::localMobilityHole(
                          mat_W_, T_C, Nf_W, E_face_W)
                    : mu_p_bulk;
                const double mu_N = (Nf_N > 1.0e10)
                    ? PhysicsEngine::localMobilityHole(
                          mat_N_, T_C, Nf_N, E_face_N)
                    : mu_p_bulk;
                const double mu_S = (Nf_S > 1.0e10)
                    ? PhysicsEngine::localMobilityHole(
                          mat_S_, T_C, Nf_S, E_face_S)
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
                const double n_i_C = static_cast<double>(m_ni_local[k]);
                const double U     = PE::netRecombination(
                    n_C, p_C, n_i_C, mat_C);
                const double Em    = electricFieldMagAt(i, j);
                const double G_bt  = PE::kaneBTBT(
                    Em, mat_C.A_kane, mat_C.B_kane, mat_C.btbt_isDirect);
                const double RmG   = U - G_bt;
                (void)mat; (void)T_lattice;

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

    // Phase 5: refresh per-cell n_i / chi-shift / eps_r cache once per
    // Gummel call.  Inner solvers read from these float buffers.
    refreshLocalParamCache(V_T);

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
// Wachutka local heat equation (FTCS, sub-stepped)            [Phase 5]
// -----------------------------------------------------------------------------
// Solves the lattice-temperature equation
//
//     rho Cp dT/dt = div(kappa grad T) + H(x,y)
//
// on the painter grid.  Heat source H(x,y) decomposes into:
//
//   H_J = sigma(x) * |E(x)|^2                     [Joule]
//       = q (n mu_n + p mu_p) |E|^2
//   H_R = (R_SRH + R_Aug - G_BTBT) * (E_g + 3 k_B T) * q   [Recombination]
//
// Discretization: cell-centred finite volume with harmonic-mean face
// thermal conductivities so the heterogeneous-medium flux is rigorous
// (kappa varies across material brushes):
//
//   bar_kappa_{i+1/2} = 2 kappa_C kappa_{C+1} / (kappa_C + kappa_{C+1})
//
// Time stepping: forward Euler (FTCS) with the Fourier-number CFL bound
//
//   alpha_eff = max_x (kappa_face_sum / rho_cp(x))
//   dt_sub <= 0.24 * h^2 / alpha_eff
//
// Sub-stepping covers the full requested wall-clock interval; if more
// than 256 sub-steps would be required we cap and let the heat lag the
// carriers by a frame (visually obvious; physically harmless).
//
// All buffers (m_T, m_T_next, m_T_old, m_H_gen, profile lookups via
// m_mat) are pre-allocated -- strict zero-allocation hot path.
//
// References: Wachutka, IEEE Trans. CAD 9 (1990) 1141; Selberherr Sec. 4.5;
//             Sze Sec. 6.6; Lindefelt, J. Appl. Phys. 75 (1994) 942.
// =============================================================================
float DriftDiffusion::solveHeatEquation(
    double n_i_ref, double V_T,
    double mu_n_ref, double mu_p_ref,
    double dt_seconds) noexcept
{
    (void)n_i_ref;
    if (dt_seconds <= 0.0) return 0.0f;
    if (V_T <= 0.0)        return 0.0f;

    // Cache refresh first: H_R reads m_ni_local; we must use the latest
    // local-T n_i, otherwise the Wachutka feedback lags one step.
    refreshLocalParamCache(V_T);

    constexpr double q  = 1.602176634e-19;     // [C]
    constexpr double kB = 8.617333262e-5;      // [eV/K]
    const double h  = static_cast<double>(m_cell_pitch_cm);
    const double h2 = h * h;

    // -------------------------------------------------------------------
    // 1) Compute H(x,y) once per call and freeze it for the sub-stepping.
    //    We do not iterate H/T self-consistently inside this call (that's
    //    what the outer Gummel + heat loop does in the main loop).
    // -------------------------------------------------------------------
    for (int j = 0; j < m_H; ++j) {
        for (int i = 0; i < m_W; ++i) {
            const std::size_t k = idx(i, j);
            const auto& mat_C  = profileAt(i, j);
            const double T_C   = static_cast<double>(m_T[k]);

            // Joule: H_J = sigma * |E|^2  with local sigma at the cell.
            const double n_C   = static_cast<double>(m_n_dens[k]);
            const double p_C   = static_cast<double>(m_p_dens[k]);
            const double E_mag = electricFieldMagAt(i, j);

            const double N_tot = static_cast<double>(m_Nd[k] + m_Na[k]);
            const double mu_n_C = (N_tot > 1.0e10)
                ? PhysicsEngine::localMobilityElectron(
                      mat_C, T_C, N_tot, E_mag)
                : mu_n_ref;
            const double mu_p_C = (N_tot > 1.0e10)
                ? PhysicsEngine::localMobilityHole(
                      mat_C, T_C, N_tot, E_mag)
                : mu_p_ref;

            const double sigma_C = q * (n_C * mu_n_C + p_C * mu_p_C);
            const double H_J     = sigma_C * E_mag * E_mag;       // W/cm^3

            // Recombination: H_R = (U_net - G_BTBT) * (E_g + 3 kT) * q
            //   units: (cm^-3 s^-1) * eV * J/eV = J / (cm^3 s) = W/cm^3
            //
            // U_net includes radiative now. In a regular diode/BJT the
            // emitted photons are reabsorbed locally (high-index Si traps
            // ~96% of internal light), so they show up as Joule-like
            // lattice heat. For an LED with engineered photon extraction
            // the radiative term should NOT be heat-coupled -- a future
            // PhotonTransport pass can subtract eta_ext * R_rad from H_R.
            const double n_i_C = static_cast<double>(m_ni_local[k]);
            const double U     = PhysicsEngine::netRecombination(
                n_C, p_C, n_i_C, mat_C);
            const double G_bt  = PhysicsEngine::kaneBTBT(
                E_mag, mat_C.A_kane, mat_C.B_kane, mat_C.btbt_isDirect);
            const double Eg_C  = PhysicsEngine::bandgapAt(mat_C, T_C);   // eV
            const double E_th  = 3.0 * kB * T_C;                          // eV
            const double H_R   = (U - G_bt) * (Eg_C + E_th) * q;

            m_H_gen[k] = static_cast<float>(H_J + H_R);
        }
    }

    // -------------------------------------------------------------------
    // 2) Compute global CFL bound.  The strict stability condition for
    //    FTCS with non-uniform kappa is
    //
    //       dt_sub * (sum of kappa_face) / (h^2 * rho_cp_C) <= 1
    //
    //    so we scan once for the worst cell.  In a uniform reference
    //    material this reduces to dt * 4 kappa / (h^2 rho_cp) <= 1, i.e.,
    //    Fourier <= 1/4.
    // -------------------------------------------------------------------
    double diff_max = 1.0e-6;     // [1/s] avoid div-by-zero
    for (int j = 0; j < m_H; ++j) {
        for (int i = 0; i < m_W; ++i) {
            const auto& mat = profileAt(i, j);
            const double rho_cp_C = std::max(mat.rho_cp, 1.0e-12);
            // Worst case: all four neighbours have the same kappa as
            // this cell (homogeneous patch). Cheap, slightly conservative.
            const double diff = 4.0 * mat.kappa / (rho_cp_C * h2);
            if (diff > diff_max) diff_max = diff;
        }
    }

    constexpr double SAFETY = 0.24;     // sub-Fourier safety margin
    const double dt_max = SAFETY / diff_max;
    int  n_sub = static_cast<int>(std::ceil(dt_seconds / dt_max));
    if (n_sub < 1)   n_sub = 1;
    if (n_sub > 256) n_sub = 256;       // cap; UI sees a one-frame lag
    const double dt_sub = dt_seconds / static_cast<double>(n_sub);

    // -------------------------------------------------------------------
    // 3) FTCS sub-stepping. Dirichlet edges at T_ambient (heat sink).
    // -------------------------------------------------------------------
    std::copy(m_T.begin(), m_T.end(), m_T_old.begin());
    float dT_max = 0.0f;

    for (int s = 0; s < n_sub; ++s) {
        for (int j = 0; j < m_H; ++j) {
            for (int i = 0; i < m_W; ++i) {
                const std::size_t k = idx(i, j);

                if (i == 0 || j == 0 || i == m_W - 1 || j == m_H - 1) {
                    m_T_next[k] = m_T_ambient;
                    continue;
                }

                const int ip = i + 1, im = i - 1;
                const int jp = j + 1, jm = j - 1;
                const std::size_t kE = idx(ip, j), kW = idx(im, j);
                const std::size_t kN = idx(i, jp), kS = idx(i, jm);

                const auto& mat_C = profileAt(i, j);
                const double rho_cp_C = std::max(mat_C.rho_cp, 1.0e-12);

                const double kappa_C = mat_C.kappa;
                const double kappa_E = profileAt(ip, j).kappa;
                const double kappa_W = profileAt(im, j).kappa;
                const double kappa_N = profileAt(i,  jp).kappa;
                const double kappa_S = profileAt(i,  jm).kappa;

                const double kf_E = 2.0 * kappa_C * kappa_E
                                  / (kappa_C + kappa_E);
                const double kf_W = 2.0 * kappa_C * kappa_W
                                  / (kappa_C + kappa_W);
                const double kf_N = 2.0 * kappa_C * kappa_N
                                  / (kappa_C + kappa_N);
                const double kf_S = 2.0 * kappa_C * kappa_S
                                  / (kappa_C + kappa_S);

                const double T_C = static_cast<double>(m_T[k]);
                const double div_term = (
                      kf_E * (m_T[kE] - T_C)
                    + kf_W * (m_T[kW] - T_C)
                    + kf_N * (m_T[kN] - T_C)
                    + kf_S * (m_T[kS] - T_C)
                ) / h2;

                const double H = static_cast<double>(m_H_gen[k]);
                double T_new = T_C + dt_sub * (div_term + H) / rho_cp_C;
                T_new = std::clamp(T_new, 50.0, 1500.0);

                const float dT = std::abs(
                    static_cast<float>(T_new - T_C));
                if (dT > dT_max) dT_max = dT;

                m_T_next[k] = static_cast<float>(T_new);
            }
        }
        m_T.swap(m_T_next);
    }

    m_dT_max_last = dT_max;
    return dT_max;
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

    // J at V + dV  (residual discarded -- we only need the converged J)
    setAppliedBias(V_save + static_cast<float>(dV));
    (void)solveGummel(n_i, V_T, epsilon_r, mu_n, mu_p, mat, T_lattice,
                      OUTER, POISS, CONT, 0.85, 1.0);
    const double J_plus  = terminalCurrentDensity(n_i, V_T, mu_n, mu_p);

    // restore mid-point
    std::copy(save_n.begin(),     save_n.end(),     m_n_dens.begin());
    std::copy(save_p.begin(),     save_p.end(),     m_p_dens.begin());
    std::copy(save_psi.begin(),   save_psi.end(),   m_psi.begin());
    std::copy(save_phi_n.begin(), save_phi_n.end(), m_phi_n.begin());
    std::copy(save_phi_p.begin(), save_phi_p.end(), m_phi_p.begin());

    // J at V - dV  (residual discarded -- we only need the converged J)
    setAppliedBias(V_save - static_cast<float>(dV));
    (void)solveGummel(n_i, V_T, epsilon_r, mu_n, mu_p, mat, T_lattice,
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


// =============================================================================
// Phase 7 -- Fully-coupled Newton-Raphson with JFNK + BiCGSTAB
// -----------------------------------------------------------------------------
// Layout of the 3N-long state vector x:
//
//     [ psi(0) ... psi(N-1) | n(0) ... n(N-1) | p(0) ... p(N-1) ]
//
//   k_psi = idx(i,j)
//   k_n   = N + idx(i,j)
//   k_p   = 2N + idx(i,j)
//
// The same layout is used for the residual F(x). Contact cells contribute
// pure Dirichlet residuals (x - x_bc); interior cells emit the standard
// drift-diffusion residuals.
//
// References: Selberherr Sec. 7.4 (coupled Newton); van der Vorst, SIAM
// J. Sci. Stat. Comput. 13 (1992) 631 (BiCGSTAB); Knoll & Keyes,
// J. Comput. Phys. 193 (2004) 357 (JFNK survey).
// =============================================================================
namespace {

// Inner-product / norm helpers. Pure-functional, no allocation.
[[nodiscard]] double dot(std::span<const double> a,
                         std::span<const double> b) noexcept {
    double s = 0.0;
    const std::size_t n = std::min(a.size(), b.size());
    for (std::size_t i = 0; i < n; ++i) s += a[i] * b[i];
    return s;
}

[[nodiscard]] double norm2(std::span<const double> a) noexcept {
    return std::sqrt(dot(a, a));
}

[[nodiscard]] double normInf(std::span<const double> a) noexcept {
    double m = 0.0;
    for (double v : a) m = std::max(m, std::abs(v));
    return m;
}

void axpy(double a, std::span<const double> x,
          std::span<double> y) noexcept {
    const std::size_t n = std::min(x.size(), y.size());
    for (std::size_t i = 0; i < n; ++i) y[i] += a * x[i];
}

} // namespace


void DriftDiffusion::packState(std::span<double> x) const noexcept {
    const std::size_t N = static_cast<std::size_t>(m_W) *
                          static_cast<std::size_t>(m_H);
    if (x.size() < 3 * N) return;
    for (std::size_t k = 0; k < N; ++k) {
        x[k]         = static_cast<double>(m_psi[k]);
        x[N + k]     = static_cast<double>(m_n_dens[k]);
        x[2 * N + k] = static_cast<double>(m_p_dens[k]);
    }
}

void DriftDiffusion::unpackState(std::span<const double> x) noexcept {
    const std::size_t N = static_cast<std::size_t>(m_W) *
                          static_cast<std::size_t>(m_H);
    if (x.size() < 3 * N) return;
    for (std::size_t k = 0; k < N; ++k) {
        m_psi   [k] = static_cast<float>(x[k]);
        // Density floor: BiCGSTAB intermediate states may dip below 0,
        // but the committed state must stay physical for the next
        // residual / Gummel call.
        m_n_dens[k] = static_cast<float>(std::max(x[N     + k], 0.0));
        m_p_dens[k] = static_cast<float>(std::max(x[2 * N + k], 0.0));
    }
}


// =============================================================================
// Residual evaluator -- F(x) for the full drift-diffusion system
// -----------------------------------------------------------------------------
// Per non-contact cell (i,j):
//
//   F_psi = (eps_s / h^2) [psi_E + psi_W + psi_N + psi_S - 4 psi_C]
//         + q (p_C - n_C + Nd - Na)                              [Sze 2.4]
//
//   F_n   = (V_T / h^2) [ mu_E B(x_E) n_E + mu_W B(-x_W) n_W
//                       + mu_N B(x_N) n_N + mu_S B(-x_S) n_S
//                       - (mu_E B(-x_E) + mu_W B(x_W)
//                        + mu_N B(-x_N) + mu_S B(x_S)) n_C ]
//         + (G_BTBT - U_net)                                     [Selberherr 6.2]
//
//   F_p   = same with B(x) <-> B(-x) (opposite drift sign)
//         + (G_BTBT - U_net)
//
// For transient (BE):  F_n -= (n_C - n_old) / dt;  same for F_p.
// At contact cells:  F = x - x_bc.
//
// Pure-functional: reads only x + m_Nd, m_Na, m_contact, m_n_dens_old,
// m_p_dens_old, m_V_bias, m_dt, m_transient_enabled, m_cell_pitch_cm.
// =============================================================================
void DriftDiffusion::computeResidual(
    std::span<const double> x,
    std::span<double>       F,
    double n_i, double V_T, double epsilon_r,
    double mu_n_bulk, double mu_p_bulk,
    const material::Profile& mat,
    double T_lattice) const noexcept
{
    const std::size_t N = static_cast<std::size_t>(m_W) *
                          static_cast<std::size_t>(m_H);
    if (x.size() < 3 * N || F.size() < 3 * N) return;
    if (V_T <= 0.0 || n_i <= 0.0) return;

    constexpr double eps_0_local = 8.854187817e-14;     // F/cm
    constexpr double q           = 1.602176634e-19;     // C
    const double eps_s = eps_0_local * epsilon_r;
    const double h     = static_cast<double>(m_cell_pitch_cm);
    const double h2    = h * h;
    const double inv_dt = (m_transient_enabled && m_dt > 0.0)
                        ? (1.0 / m_dt) : 0.0;

    using PE = PhysicsEngine;

    // Density floor used only for recombination & B-coefficient args to
    // avoid log/exp of negatives during transient JFNK perturbations.
    constexpr double DENS_FLOOR = 1.0;

    for (int j = 0; j < m_H; ++j) {
        for (int i = 0; i < m_W; ++i) {
            const std::size_t k = idx(i, j);

            // ---- Contact (Dirichlet) ---------------------------------------
            if (m_contact[k] != 0u) {
                const double V_metal = (m_contact[k] == 1u)
                    ? static_cast<double>(m_V_bias) : 0.0;
                const double Nd_C = std::max(double{m_Nd[k]}, 1.0);
                const double Na_C = std::max(double{m_Na[k]}, 1.0);
                const bool   isAnode = (m_contact[k] == 1u);
                const double psi_bc = isAnode
                    ? V_metal - V_T * std::log(Na_C / n_i)
                    : V_metal + V_T * std::log(Nd_C / n_i);
                const double n_bc = isAnode ? (n_i * n_i) / Na_C : Nd_C;
                const double p_bc = isAnode ? Na_C : (n_i * n_i) / Nd_C;

                F[k]         = x[k]         - psi_bc;
                F[N + k]     = x[N + k]     - n_bc;
                F[2 * N + k] = x[2 * N + k] - p_bc;
                continue;
            }

            // ---- Neighbour indices (Neumann edge clamp) --------------------
            const int ip = std::min(i + 1, m_W - 1);
            const int im = std::max(i - 1, 0);
            const int jp = std::min(j + 1, m_H - 1);
            const int jm = std::max(j - 1, 0);
            const std::size_t kE = idx(ip, j);
            const std::size_t kW = idx(im, j);
            const std::size_t kN = idx(i, jp);
            const std::size_t kS = idx(i, jm);

            // ---- Decode x at this stencil ----------------------------------
            const double psi_C = x[k];
            const double psi_E = x[kE];
            const double psi_W = x[kW];
            const double psi_N = x[kN];
            const double psi_S = x[kS];
            const double n_C   = std::max(x[N + k],  DENS_FLOOR);
            const double n_E   = std::max(x[N + kE], DENS_FLOOR);
            const double n_W   = std::max(x[N + kW], DENS_FLOOR);
            const double n_N   = std::max(x[N + kN], DENS_FLOOR);
            const double n_S   = std::max(x[N + kS], DENS_FLOOR);
            const double p_C   = std::max(x[2 * N + k],  DENS_FLOOR);
            const double p_E   = std::max(x[2 * N + kE], DENS_FLOOR);
            const double p_W   = std::max(x[2 * N + kW], DENS_FLOOR);
            const double p_N   = std::max(x[2 * N + kN], DENS_FLOOR);
            const double p_S   = std::max(x[2 * N + kS], DENS_FLOOR);

            // ---- F_psi: Poisson 5-point + charge density -------------------
            const double rho =
                p_C - n_C
                + static_cast<double>(m_Nd[k]) - static_cast<double>(m_Na[k]);
            F[k] = (eps_s / h2) *
                   (psi_E + psi_W + psi_N + psi_S - 4.0 * psi_C)
                 + q * rho;

            // ---- SG coefficients (one B(+x) and one B(-x) per face) --------
            const double xE = (psi_E - psi_C) / V_T;
            const double xW = (psi_C - psi_W) / V_T;
            const double xN = (psi_N - psi_C) / V_T;
            const double xS = (psi_C - psi_S) / V_T;
            const double bE  = PE::bernoulli( xE);
            const double bnE = PE::bernoulli(-xE);
            const double bW  = PE::bernoulli( xW);
            const double bnW = PE::bernoulli(-xW);
            const double bN  = PE::bernoulli( xN);
            const double bnN = PE::bernoulli(-xN);
            const double bS  = PE::bernoulli( xS);
            const double bnS = PE::bernoulli(-xS);

            // ---- Per-face mobility (Matthiessen + Caughey-Thomas) ----------
            // |E|_face from one-sided psi difference.  Local doping at face
            // is the arithmetic mean of the two adjoining cells.
            const double E_E = std::abs(psi_E - psi_C) / h;
            const double E_W = std::abs(psi_C - psi_W) / h;
            const double E_N = std::abs(psi_N - psi_C) / h;
            const double E_S = std::abs(psi_C - psi_S) / h;
            const double Nf_E = 0.5 * (m_Nd[k] + m_Nd[kE]
                                      + m_Na[k] + m_Na[kE]);
            const double Nf_W = 0.5 * (m_Nd[k] + m_Nd[kW]
                                      + m_Na[k] + m_Na[kW]);
            const double Nf_N = 0.5 * (m_Nd[k] + m_Nd[kN]
                                      + m_Na[k] + m_Na[kN]);
            const double Nf_S = 0.5 * (m_Nd[k] + m_Nd[kS]
                                      + m_Na[k] + m_Na[kS]);
            const double muE_n = (Nf_E > 1.0e10)
                ? PE::localMobilityElectron(mat, T_lattice, Nf_E, E_E)
                : mu_n_bulk;
            const double muW_n = (Nf_W > 1.0e10)
                ? PE::localMobilityElectron(mat, T_lattice, Nf_W, E_W)
                : mu_n_bulk;
            const double muN_n = (Nf_N > 1.0e10)
                ? PE::localMobilityElectron(mat, T_lattice, Nf_N, E_N)
                : mu_n_bulk;
            const double muS_n = (Nf_S > 1.0e10)
                ? PE::localMobilityElectron(mat, T_lattice, Nf_S, E_S)
                : mu_n_bulk;
            const double muE_p = (Nf_E > 1.0e10)
                ? PE::localMobilityHole(mat, T_lattice, Nf_E, E_E)
                : mu_p_bulk;
            const double muW_p = (Nf_W > 1.0e10)
                ? PE::localMobilityHole(mat, T_lattice, Nf_W, E_W)
                : mu_p_bulk;
            const double muN_p = (Nf_N > 1.0e10)
                ? PE::localMobilityHole(mat, T_lattice, Nf_N, E_N)
                : mu_p_bulk;
            const double muS_p = (Nf_S > 1.0e10)
                ? PE::localMobilityHole(mat, T_lattice, Nf_S, E_S)
                : mu_p_bulk;

            // ---- F_n: electron continuity in residual form -----------------
            const double lap_n = (V_T / h2) * (
                muE_n * bE  * n_E + muW_n * bnW * n_W
              + muN_n * bN  * n_N + muS_n * bnS * n_S
              - (muE_n * bnE + muW_n * bW
               + muN_n * bnN + muS_n * bS) * n_C
            );

            // Spatially varying local n_i (Phase 5) is supported via
            // m_ni_local when present; fall back to the global n_i if not
            // currently meaningful (uniform material).
            const double n_i_C = (k < m_ni_local.size())
                ? static_cast<double>(m_ni_local[k]) : n_i;
            const double U  = PE::netRecombination(n_C, p_C, n_i_C, mat);
            const double Em = std::sqrt(
                std::pow((psi_E - psi_W) / (2.0 * h), 2) +
                std::pow((psi_N - psi_S) / (2.0 * h), 2));
            const double G_bt = PE::kaneBTBT(
                Em, mat.A_kane, mat.B_kane, mat.btbt_isDirect);

            double F_n = lap_n + (G_bt - U);
            if (inv_dt > 0.0) {
                F_n -= (n_C - static_cast<double>(m_n_dens_old[k])) * inv_dt;
            }
            F[N + k] = F_n;

            // ---- F_p: hole continuity (B(x) <-> B(-x) swap) ---------------
            const double lap_p = (V_T / h2) * (
                muE_p * bnE * p_E + muW_p * bW  * p_W
              + muN_p * bnN * p_N + muS_p * bS  * p_S
              - (muE_p * bE  + muW_p * bnW
               + muN_p * bN  + muS_p * bnS) * p_C
            );
            double F_p = lap_p + (G_bt - U);
            if (inv_dt > 0.0) {
                F_p -= (p_C - static_cast<double>(m_p_dens_old[k])) * inv_dt;
            }
            F[2 * N + k] = F_p;
        }
    }
}


// =============================================================================
// JFNK matrix-vector product:  Jv ~ (F(x + eps v) - F(x)) / eps
// -----------------------------------------------------------------------------
// eps_FD = sqrt(eps_machine) * (1 + ||x||_2) / ||v||_2   (Knoll/Keyes 4.4)
// chosen to balance roundoff against truncation error.
// =============================================================================
void DriftDiffusion::jfnkMatvec(
    std::span<const double> v,
    std::span<double>       Jv,
    std::span<const double> x_base,
    std::span<const double> F_base,
    double n_i, double V_T, double epsilon_r,
    double mu_n_bulk, double mu_p_bulk,
    const material::Profile& mat,
    double T_lattice) noexcept
{
    const std::size_t n3 = x_base.size();
    if (v.size() < n3 || Jv.size() < n3
        || m_nk_x_pert.size() < n3 || m_nk_F_pert.size() < n3
        || F_base.size() < n3) return;

    const double v_norm = norm2(v);
    if (v_norm <= 0.0) {
        for (std::size_t i = 0; i < n3; ++i) Jv[i] = 0.0;
        return;
    }
    const double x_norm = norm2(x_base);
    const double eps_mach = 1.490116e-8;            // sqrt(DBL_EPSILON)
    const double eps_fd   = eps_mach * (1.0 + x_norm) / v_norm;

    for (std::size_t i = 0; i < n3; ++i)
        m_nk_x_pert[i] = x_base[i] + eps_fd * v[i];

    computeResidual(
        std::span<const double>(m_nk_x_pert.data(), n3),
        std::span<double>      (m_nk_F_pert.data(), n3),
        n_i, V_T, epsilon_r,
        mu_n_bulk, mu_p_bulk,
        mat, T_lattice);

    const double inv_eps = 1.0 / eps_fd;
    for (std::size_t i = 0; i < n3; ++i)
        Jv[i] = (m_nk_F_pert[i] - F_base[i]) * inv_eps;
}


// =============================================================================
// BiCGSTAB linear solver (van der Vorst 1992)
// -----------------------------------------------------------------------------
// Solves  J dx = b  (J accessed only via JFNK matvec).
// Standard pseudocode, with breakdown guards (rho ~ 0, omega ~ 0).
// Initial guess x_out = 0 -- safe and avoids state coupling between Newton
// steps.  Returns the iteration count at exit.
// =============================================================================
int DriftDiffusion::bicgstabSolve(
    std::span<const double> b,
    std::span<double>       x_out,
    int max_iter, double tol,
    std::span<const double> x_base,
    std::span<const double> F_base,
    double n_i, double V_T, double epsilon_r,
    double mu_n_bulk, double mu_p_bulk,
    const material::Profile& mat,
    double T_lattice) noexcept
{
    const std::size_t n3 = b.size();
    if (x_out.size() < n3) return 0;

    auto R  = std::span<double>(m_bicg_r .data(), n3);
    auto Rh = std::span<double>(m_bicg_rh.data(), n3);
    auto P  = std::span<double>(m_bicg_p .data(), n3);
    auto V  = std::span<double>(m_bicg_v .data(), n3);
    auto S  = std::span<double>(m_bicg_s .data(), n3);
    auto T  = std::span<double>(m_bicg_t .data(), n3);

    // x_out = 0; r = b - A*x_out = b
    for (std::size_t i = 0; i < n3; ++i) {
        x_out[i] = 0.0;
        R [i]    = b[i];
        Rh[i]    = b[i];          // shadow residual = initial residual
        P [i]    = 0.0;
        V [i]    = 0.0;
    }

    const double b_norm = norm2(b);
    if (b_norm <= 0.0) return 0;
    const double abs_tol = tol * b_norm;

    double rho_old = 1.0, alpha = 1.0, omega_bcg = 1.0;
    int iter = 0;

    for (iter = 1; iter <= max_iter; ++iter) {
        const double rho = dot(Rh, R);
        if (std::abs(rho) < 1.0e-30) break;   // breakdown

        const double beta = (rho / rho_old) * (alpha / omega_bcg);
        // p = r + beta * (p - omega * v)
        for (std::size_t i = 0; i < n3; ++i)
            P[i] = R[i] + beta * (P[i] - omega_bcg * V[i]);

        jfnkMatvec(P, V,
                   x_base, F_base,
                   n_i, V_T, epsilon_r,
                   mu_n_bulk, mu_p_bulk, mat, T_lattice);

        const double rh_v = dot(Rh, V);
        if (std::abs(rh_v) < 1.0e-30) break;
        alpha = rho / rh_v;

        // s = r - alpha * v
        for (std::size_t i = 0; i < n3; ++i)
            S[i] = R[i] - alpha * V[i];

        if (norm2(S) < abs_tol) {
            // Early convergence: x += alpha * p; done.
            axpy(alpha, P, x_out);
            ++iter;
            break;
        }

        jfnkMatvec(S, T,
                   x_base, F_base,
                   n_i, V_T, epsilon_r,
                   mu_n_bulk, mu_p_bulk, mat, T_lattice);

        const double t_t = dot(T, T);
        if (t_t < 1.0e-30) break;
        omega_bcg = dot(T, S) / t_t;

        // x += alpha * p + omega * s
        for (std::size_t i = 0; i < n3; ++i)
            x_out[i] += alpha * P[i] + omega_bcg * S[i];

        // r = s - omega * t
        for (std::size_t i = 0; i < n3; ++i)
            R[i] = S[i] - omega_bcg * T[i];

        if (norm2(R) < abs_tol) break;
        if (std::abs(omega_bcg) < 1.0e-30) break;
        rho_old = rho;
    }
    return iter;
}


// =============================================================================
// solveNewton -- public outer loop
// -----------------------------------------------------------------------------
// 1. Seed densities (Boltzmann) on first call.
// 2. Apply contact Dirichlet BCs to (psi, n, p).
// 3. Pack state -> x.
// 4. Newton loop:
//    a. F = compute_residual(x)
//    b. If ||F||_inf < newton_tol -> done.
//    c. Solve J dx = -F via BiCGSTAB + JFNK.
//    d. Damped + psi-clamped update: x += alpha * dx, |dpsi| <= V_T_clamp.
// 5. Unpack x -> state, refresh quasi-Fermi for the BandView.
//
// Returns the final ||F||_inf at convergence (or after newton_iters
// exhausted).  Designed so the call site can swap solveGummel <->
// solveNewton freely; both produce the same (psi, n, p, phi_n, phi_p)
// observable state.
// =============================================================================
double DriftDiffusion::solveNewton(
    double n_i, double V_T, double epsilon_r,
    double mu_n_bulk, double mu_p_bulk,
    const material::Profile& mat,
    double T_lattice,
    int    newton_iters,
    int    bicg_iters,
    double bicg_tol,
    double newton_tol,
    double damping) noexcept
{
    const std::size_t N  = static_cast<std::size_t>(m_W) *
                           static_cast<std::size_t>(m_H);
    const std::size_t n3 = 3 * N;
    if (n3 == 0 || m_nk_x.size() < n3) return 0.0;

    if (!m_dens_seeded) {
        applyContactBoundaries(n_i, V_T);
        seedDensitiesFromBoltzmann(n_i, V_T);
    }
    applyContactBoundaries(n_i, V_T);

    // ---- Pack initial state -----------------------------------------------
    auto X  = std::span<double>(m_nk_x .data(), n3);
    auto DX = std::span<double>(m_nk_dx.data(), n3);
    auto F  = std::span<double>(m_nk_F .data(), n3);
    packState(X);

    // Psi-clamp: never accept a |Delta psi| > 1 V per Newton step.  Keeps
    // the Boltzmann exponents (~psi/V_T) from drifting into the regime
    // where the linearisation is meaningless.
    constexpr double PSI_CLAMP = 1.0;

    double F_inf = 0.0;
    int    last_bicg = 0;
    int    k_done    = 0;

    for (int k = 0; k < newton_iters; ++k) {
        // ---- F = F(x) -----------------------------------------------------
        computeResidual(X, F,
                        n_i, V_T, epsilon_r,
                        mu_n_bulk, mu_p_bulk,
                        mat, T_lattice);
        F_inf = normInf(F);
        k_done = k + 1;
        if (F_inf < newton_tol) break;

        // ---- b = -F -------------------------------------------------------
        // We pass -F as the RHS to BiCGSTAB; reuse F_pert as scratch.
        auto NF = std::span<double>(m_nk_F_pert.data(), n3);
        for (std::size_t i = 0; i < n3; ++i) NF[i] = -F[i];

        // ---- Solve J dx = -F ---------------------------------------------
        last_bicg = bicgstabSolve(
            std::span<const double>(NF.data(), n3),
            DX, bicg_iters, bicg_tol,
            std::span<const double>(X.data(), n3),
            std::span<const double>(F.data(), n3),
            n_i, V_T, epsilon_r,
            mu_n_bulk, mu_p_bulk, mat, T_lattice);

        // ---- Psi-clamp the update ----------------------------------------
        // First N entries of dx are dpsi; saturate each component before
        // applying.  Damping is applied uniformly afterwards.
        for (std::size_t i = 0; i < N; ++i) {
            if      (DX[i] >  PSI_CLAMP) DX[i] =  PSI_CLAMP;
            else if (DX[i] < -PSI_CLAMP) DX[i] = -PSI_CLAMP;
        }

        // ---- Apply damped update -----------------------------------------
        for (std::size_t i = 0; i < n3; ++i)
            X[i] += damping * DX[i];

        // Re-apply contacts so Dirichlet rows stay exact.
        unpackState(std::span<const double>(X.data(), n3));
        applyContactBoundaries(n_i, V_T);
        packState(X);
    }

    // ---- Commit to live buffers + refresh quasi-Fermi ---------------------
    unpackState(std::span<const double>(X.data(), n3));
    refreshQuasiFermiFromDensities(n_i, V_T);

    // ---- Bookkeeping for the UI -------------------------------------------
    m_nk_last_resid   = F_inf;
    m_nk_last_iters   = k_done;
    m_bicg_last_iters = last_bicg;
    return F_inf;
}
