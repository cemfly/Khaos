#include "DriftDiffusion.hpp"

#include <algorithm>
#include <cmath>


// =============================================================================
// Construction
// =============================================================================
DriftDiffusion::DriftDiffusion(int gridW, int gridH)
    : m_W(gridW), m_H(gridH)
    , m_n     (static_cast<std::size_t>(gridW * gridH), 0.0f)
    , m_n_next(static_cast<std::size_t>(gridW * gridH), 0.0f)
    , m_G     (static_cast<std::size_t>(gridW * gridH), 0.0f)
{}


// =============================================================================
// Material-dependent integrator tuning.
//
//   D_scale : relative diffusion coefficient (Si = 1.0 baseline). GaAs ~ 6
//             due to higher electron mobility; Ge ~ 3.
//   tau_seconds : minority-carrier lifetime (visualisation units; we keep
//                 it on a scale of seconds so spots fade over a few seconds
//                 at 60 fps).
// =============================================================================
void DriftDiffusion::configureForMaterial(double D_scale, double tau_seconds) {
    // Map relative D into a stable FTCS coefficient.
    //   For an unbounded 2D explicit scheme with grid spacing h = 1 cell,
    //   stability requires  D * dt < 0.25.  We pick D = 0.04 * D_scale and
    //   integrate at variable dt, capping the per-step coefficient below
    //   that limit inside step().
    m_D   = static_cast<float>(0.04 * D_scale);
    m_tau = static_cast<float>(std::max(0.05, tau_seconds));
}


// =============================================================================
// Source placement
// =============================================================================
void DriftDiffusion::addSource(float u, float v, float intensity, float sigma) {
    u = std::clamp(u, 0.0f, 1.0f);
    v = std::clamp(v, 0.0f, 1.0f);
    const float i0 = u * (m_W - 1);
    const float j0 = v * (m_H - 1);
    const float two_sigma_sq_inv = 1.0f / (2.0f * sigma * sigma * (m_W * m_H));

    for (int j = 0; j < m_H; ++j) {
        for (int i = 0; i < m_W; ++i) {
            const float du = (i - i0);
            const float dv = (j - j0);
            const float r2 = du * du + dv * dv;
            const float weight = std::exp(-r2 * two_sigma_sq_inv);
            m_G[idx(i, j)] += intensity * weight;
        }
    }
}

void DriftDiffusion::clear() {
    std::fill(m_n.begin(),      m_n.end(),      0.0f);
    std::fill(m_n_next.begin(), m_n_next.end(), 0.0f);
    std::fill(m_G.begin(),      m_G.end(),      0.0f);
}


// =============================================================================
// FTCS step
//
//   n_next[i,j] = n[i,j]
//               + D_eff * (n[i+1,j] + n[i-1,j] + n[i,j+1] + n[i,j-1]
//                          - 4 n[i,j])
//               + G[i,j] * dt
//               - n[i,j] * dt / tau
//
//   D_eff is clamped to <= 0.24 per step to retain unconditional stability
//   regardless of frame-rate variations.  Reflective Neumann boundaries
//   (carriers cannot leave the visualization).
// =============================================================================
void DriftDiffusion::step(float dt_seconds) {
    const float D_eff = std::min(0.24f, m_D * dt_seconds * 60.0f); // tuned per 60fps
    const float decay = dt_seconds / m_tau;

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

            float v = c + D_eff * lap + m_G[idx(i, j)] * dt_seconds - c * decay;
            v = std::max(0.0f, v);
            m_n_next[idx(i, j)] = v;
        }
    }
    m_n.swap(m_n_next);
}


// =============================================================================
// Read-out helpers
// =============================================================================
float DriftDiffusion::at(int i, int j) const noexcept {
    if (i < 0 || i >= m_W || j < 0 || j >= m_H) return 0.0f;
    return m_n[idx(i, j)];
}

float DriftDiffusion::maxValue() const noexcept {
    if (m_n.empty()) return 0.0f;
    return *std::max_element(m_n.begin(), m_n.end());
}

float DriftDiffusion::meanValue() const noexcept {
    if (m_n.empty()) return 0.0f;
    double sum = 0.0;
    for (float v : m_n) sum += v;
    return static_cast<float>(sum / m_n.size());
}

double DriftDiffusion::globalExcess() const noexcept {
    return static_cast<double>(meanValue()) * static_cast<double>(m_dN_scale);
}
