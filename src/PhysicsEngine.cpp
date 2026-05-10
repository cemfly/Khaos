// =============================================================================
// PhysicsEngine.cpp
//
//   Author : dex / cemfly-april2026
//   License: MIT
// =============================================================================

#include "PhysicsEngine.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <sstream>


// =============================================================================
// Construction
// =============================================================================
PhysicsEngine::PhysicsEngine(material::Kind kind)
    : m_materialKind(kind)
    , m_material(&material::byKind(kind))
    , m_T(300.0)
    , m_dopingType(DopingType::Intrinsic)
    , m_N(0.0)
    , m_incompleteIonization(false)
    , m_mobilityModel(MobilityModel::Matthiessen)
    , m_lambda_nm(1500.0)
    , m_opticalEnabled(false)
    , m_B(0.0)
    , m_Eg(0.0), m_Nc(0.0), m_Nv(0.0), m_ni(0.0)
    , m_n(0.0),  m_p(0.0),  m_Ef(0.0)
    , m_NdPlus(0.0), m_NaMinus(0.0)
    , m_mu_n(0.0),   m_mu_p(0.0),   m_sigma(0.0)
    , m_Ephoton(0.0), m_Gopt(0.0),  m_deltaN(0.0)
    , m_R_H(0.0)
{
    recompute();
}


// =============================================================================
// Setters
// =============================================================================
void PhysicsEngine::setMaterial(material::Kind kind) {
    m_materialKind = kind;
    m_material     = &material::byKind(kind);
    recompute();
}

void PhysicsEngine::setTemperature(double T)            { m_T = std::max(1.0, T);          recompute(); }
void PhysicsEngine::setDopingType(DopingType t)         { m_dopingType = t;                 recompute(); }
void PhysicsEngine::setDopingConcentration(double N)    { m_N = std::max(0.0, N);           recompute(); }
void PhysicsEngine::setIncompleteIonization(bool e)     { m_incompleteIonization = e;       recompute(); }
void PhysicsEngine::setMobilityModel(MobilityModel m)   { m_mobilityModel = m;              recompute(); }

void PhysicsEngine::setWavelengthNm(double lambda) {
    m_lambda_nm = std::clamp(lambda, phys::lambda_min_nm, phys::lambda_max_nm);
    recompute();
}
void PhysicsEngine::setOpticalEnabled(bool e)           { m_opticalEnabled = e;             recompute(); }

void PhysicsEngine::setMagneticField(double B) {
    m_B = std::clamp(B, phys::B_min_T, phys::B_max_T);
    recompute();
}

void PhysicsEngine::setDriftDiffusionExcess(double extra) {
    m_extraDriftDN = std::max(0.0, extra);
    // No full recompute -- this only affects transport / Hall via getTotal*().
    computeTransport();
    computeHall();
}


// =============================================================================
// Auxiliary band / Fermi-Dirac helpers
// =============================================================================
double PhysicsEngine::getDonorLevel() const noexcept {
    return m_Eg - m_material->E_donor_offset;
}
double PhysicsEngine::getAcceptorLevel() const noexcept {
    return m_material->E_acceptor_offset;
}

double PhysicsEngine::fermiDirac(double E) const noexcept {
    const double kT = phys::k_B * m_T;
    const double x  = (E - m_Ef) / kT;
    if (x >  50.0) return 0.0;
    if (x < -50.0) return 1.0;
    return 1.0 / (1.0 + std::exp(x));
}

double PhysicsEngine::getIonizationFraction() const noexcept {
    if (m_dopingType == DopingType::Intrinsic || m_N <= 0.0) return 1.0;
    if (m_dopingType == DopingType::NType) return m_NdPlus  / m_N;
    return                                        m_NaMinus / m_N;
}

double PhysicsEngine::getResistivity() const noexcept {
    return (m_sigma > 0.0) ? 1.0 / m_sigma : std::numeric_limits<double>::infinity();
}

bool PhysicsEngine::isOpticallyPumped() const noexcept {
    return m_opticalEnabled && (m_Ephoton > m_Eg);
}

double PhysicsEngine::photonEnergyEv(double lambda_nm) noexcept {
    return (lambda_nm > 0.0) ? phys::hc_eVnm / lambda_nm : 0.0;
}


// =============================================================================
// Static T-evaluators (used by per-cell electrothermal feedback and tests).
//
//   bandgapAt(mat, T)         = E_g(T)  via Varshni for `mat`           [eV]
//   intrinsicCarrierAt        = n_i = sqrt(Nc Nv) exp(-Eg/2kT)          [cm^-3]
//   log10IntrinsicCarrierAt   = log10 n_i  (numerically stable form)
//
// Why a log-domain n_i evaluator?  At T = 10 K with Si the exponent
// Eg/(2kT) ~ 678, so exp(-678) underflows to *exactly zero* in double
// precision -- which makes downstream divisions blow up.  Computing
// log10(n_i) directly side-steps the underflow:
//
//     log10(n_i) = 0.5 * log10(Nc Nv) - 0.4342945 * Eg / (2 kT)
//
// Callers that only need the magnitude (graph axes, freeze-out
// thresholds, ...) can use this instead of intrinsicCarrierAt() and
// avoid silent zeros.
// =============================================================================
double PhysicsEngine::bandgapAt(const material::Profile& mat, double T) noexcept {
    if (T < 1.0) T = 1.0;
    return mat.Eg0 - (mat.varshni_a * T * T) / (T + mat.varshni_b);
}

double PhysicsEngine::log10IntrinsicCarrierAt(const material::Profile& mat,
                                              double T) noexcept
{
    if (T < 1.0) T = 1.0;
    const double Eg    = bandgapAt(mat, T);
    const double T_3_2 = std::pow(T / 300.0, 1.5);
    const double Nc    = mat.Nc_300K * T_3_2;
    const double Nv    = mat.Nv_300K * T_3_2;
    const double kT    = phys::k_B * T;
    constexpr double kLn10Inv = 0.43429448190325176;     // 1 / ln(10)
    return 0.5 * std::log10(Nc * Nv) - kLn10Inv * (Eg / (2.0 * kT));
}

double PhysicsEngine::intrinsicCarrierAt(const material::Profile& mat,
                                         double T) noexcept
{
    const double l10 = log10IntrinsicCarrierAt(mat, T);
    // IEEE754 double range: [~1e-308, ~1e+308].  Below ~1e-300 the
    // result is unrepresentable; pin to zero rather than letting a
    // denormalised tail leak into the rest of the pipeline.
    if (l10 < -300.0) return 0.0;
    if (l10 >  300.0) return std::numeric_limits<double>::infinity();
    return std::pow(10.0, l10);
}


// =============================================================================
// Caughey-Thomas high-field saturation
//
//   v(E) = mu_low * E / [1 + (mu_low * E / v_sat)^beta]^(1/beta)
//   mu_eff(E) = v(E) / E
//
// Reference: Caughey & Thomas, Proc. IEEE 55 (1967), tabulated in Sze
// "Physics of Semiconductor Devices" Sec. 1.5.4.
//
// At E -> 0: mu_eff -> mu_low
// At E -> inf: v -> v_sat (the limiting drift speed at high fields)
// =============================================================================
double PhysicsEngine::highFieldMobility(double mu_low,
                                        double E_field_V_per_cm,
                                        double v_sat_cm_per_s,
                                        double beta) noexcept
{
    if (E_field_V_per_cm <= 0.0) return mu_low;
    if (v_sat_cm_per_s   <= 0.0) return mu_low;

    const double r   = (mu_low * E_field_V_per_cm) / v_sat_cm_per_s;
    const double rb  = std::pow(r, beta);
    const double den = std::pow(1.0 + rb, 1.0 / beta);
    return mu_low / den;     // = v(E)/E
}


// =============================================================================
// Shockley-Read-Hall recombination through midgap traps.
//
//   R = (n p - n_i^2) / [tau_p (n + n_t) + tau_n (p + p_t)]
//
// We assume midgap traps  (E_t == E_i)  =>  n_t = p_t = n_i, so
//
//   R = (n p - n_i^2) / [tau_p (n + n_i) + tau_n (p + n_i)]
//
// Note: in equilibrium  n p = n_i^2  =>  R = 0,  satisfying detailed
// balance automatically.  The simple "n / tau" form used elsewhere
// violates this in equilibrium.
//
// Reference: Pierret, "Semiconductor Device Fundamentals" Ch. 5;
// Sze Sec. 1.5.5.
// =============================================================================
double PhysicsEngine::recombSRH(double n, double p, double n_i,
                                double tau_n, double tau_p) noexcept
{
    if (tau_n <= 0.0 || tau_p <= 0.0) return 0.0;
    const double num = n * p - n_i * n_i;
    const double den = tau_p * (n + n_i) + tau_n * (p + n_i);
    if (den <= 0.0) return 0.0;
    return num / den;
}


// =============================================================================
// Bandgap narrowing (Slotboom & de Graaff, Solid-State Electronics 19 (1976))
//
//   dE_g(N) = E_BGN * [ ln(N / N_ref) + sqrt(ln^2(N / N_ref) + 0.5) ]
//
// with E_BGN = 9 meV and N_ref = 1e17 cm^-3 -- the canonical fit for Si.
//
// For N below N_ref the function is essentially zero; for N >= 1e18 it
// becomes significant and depresses the *effective* intrinsic carrier
// to n_ie = n_i * exp(dE_g / 2kT). This matters for both Fermi level
// placement and BJT injection efficiency.
// =============================================================================
double PhysicsEngine::bandgapNarrowing(double N_per_cm3) noexcept {
    if (N_per_cm3 <= 0.0) return 0.0;
    constexpr double E_BGN = 9.0e-3;        // eV
    constexpr double N_ref = 1.0e17;        // cm^-3
    const double x = std::log(N_per_cm3 / N_ref);
    const double term = x + std::sqrt(x * x + 0.5);
    return std::max(0.0, E_BGN * term);
}

double PhysicsEngine::effectiveIntrinsicCarrier(const material::Profile& mat,
                                                double T,
                                                double N_per_cm3) noexcept
{
    if (T < 1.0) T = 1.0;
    const double ni     = intrinsicCarrierAt(mat, T);
    const double dEg    = bandgapNarrowing(N_per_cm3);
    const double kT     = phys::k_B * T;
    return ni * std::exp(dEg / (2.0 * kT));
}


// =============================================================================
// Optical absorption (Tauc-style for direct vs indirect bandgap).
//
//   alpha(hv) = A * (hv - E_g)^p / hv      [cm^-1]
//
// with p = 1/2 for direct allowed transitions (GaAs) and p = 2 for
// indirect transitions (Si, Ge). The 1/hv prefactor follows from the
// Joint Density of States derivation (Pankove Ch. 3).
//
// The amplitude A is calibrated so that mid-visible (hv ~ 2.5 eV) gives
// a textbook absorption coefficient (~ 1e4 /cm for Si, ~ 1e5 /cm for GaAs)
// without requiring per-material lookup tables.
// =============================================================================
double PhysicsEngine::absorptionCoefficient(const material::Profile& mat,
                                            double photon_eV) noexcept
{
    const double Eg = bandgapAt(mat, 300.0);
    if (photon_eV <= Eg) return 0.0;

    const double excess = photon_eV - Eg;
    const double p      = mat.isDirectBandgap ? 0.5 : 2.0;
    // Calibrated against textbook tables (Sze App. C, Pankove Ch. 3):
    //   Si  at 2.5 eV  ->  alpha ~  1e4 /cm
    //   GaAs at 2.5 eV ->  alpha ~  4e4 /cm
    // The direct-gap material is several times more absorbing in the
    // visible, even though its bandgap is wider than Si's.
    const double A = mat.isDirectBandgap ? 1.0e5 : 1.0e4;
    return A * std::pow(excess, p) / photon_eV;
}

double PhysicsEngine::penetrationDepthCm(const material::Profile& mat,
                                         double photon_eV) noexcept
{
    const double a = absorptionCoefficient(mat, photon_eV);
    if (a <= 0.0) return std::numeric_limits<double>::infinity();
    return 1.0 / a;
}


// =============================================================================
// BJT (NPN, one-sided emitter junction).
//
// Emitter injection efficiency (Sze Sec. 5.2):
//
//   gamma = 1 / (1 + (D_p N_B W_E) / (D_n N_E W_B))
//
// For a well-designed NPN, N_E >> N_B  =>  gamma -> 1.
//
// Common-emitter current gain:
//
//   beta = gamma * alpha_T / (1 - gamma * alpha_T)
//
// Early effect (base-width modulation):
//
//   I_C(V_CE) = I_C0 * (1 + V_CE / V_A)
// =============================================================================
double PhysicsEngine::bjtEmitterEfficiency(double D_n, double N_E, double W_E,
                                           double D_p, double N_B, double W_B) noexcept
{
    if (D_n <= 0.0 || D_p <= 0.0) return 0.0;
    if (N_E <= 0.0 || N_B <= 0.0) return 0.0;
    if (W_E <= 0.0 || W_B <= 0.0) return 0.0;
    const double ratio = (D_p * N_B * W_E) / (D_n * N_E * W_B);
    return 1.0 / (1.0 + ratio);
}

double PhysicsEngine::bjtCurrentGain(double gamma, double alpha_T) noexcept {
    const double a = std::clamp(gamma * alpha_T, 0.0, 0.99999);
    return a / (1.0 - a);
}

double PhysicsEngine::earlyEffectFactor(double V_CE, double V_Early) noexcept {
    if (V_Early <= 0.0) return 1.0;
    return 1.0 + V_CE / V_Early;
}


// =============================================================================
// Matthiessen mobility (parameterized by Material)
// =============================================================================
double PhysicsEngine::matthiessenMobilityElectron(
    const material::Profile& mat, double T, double N) noexcept
{
    const double Tr   = T / 300.0;
    const double mu_L = mat.mu_L_n_300 * std::pow(Tr, -1.5);
    if (N <= 0.0) return mu_L;
    const double mu_I = mat.mu_I_n_300 * std::pow(Tr, 1.5)
                                       * (mat.N_ref_matt / N);
    return 1.0 / (1.0 / mu_L + 1.0 / mu_I);
}

double PhysicsEngine::matthiessenMobilityHole(
    const material::Profile& mat, double T, double N) noexcept
{
    const double Tr   = T / 300.0;
    const double mu_L = mat.mu_L_p_300 * std::pow(Tr, -1.5);
    if (N <= 0.0) return mu_L;
    const double mu_I = mat.mu_I_p_300 * std::pow(Tr, 1.5)
                                       * (mat.N_ref_matt / N);
    return 1.0 / (1.0 / mu_L + 1.0 / mu_I);
}


// =============================================================================
// Arora model -- Si-tuned reference (legacy, for comparison plots).
// =============================================================================
namespace {
    constexpr double mu_n_min_300 = 88.0;
    constexpr double mu_n_max_300 = 1414.0;
    constexpr double Nref_n_300   = 1.26e17;
    constexpr double alpha_n_300  = 0.88;
    constexpr double mu_p_min_300 = 54.3;
    constexpr double mu_p_max_300 = 470.5;
    constexpr double Nref_p_300   = 2.35e17;
    constexpr double alpha_p_300  = 0.88;
}

double PhysicsEngine::aroraMobilityElectron(double T, double N) noexcept {
    const double Tr     = T / 300.0;
    const double mu_min = mu_n_min_300 * std::pow(Tr, -0.57);
    const double mu_max = mu_n_max_300 * std::pow(Tr, -2.33);
    const double N_ref  = Nref_n_300   * std::pow(Tr,  2.4);
    const double alpha  = alpha_n_300  * std::pow(Tr, -0.146);
    const double scale  = (N > 0.0) ? std::pow(N / N_ref, alpha) : 0.0;
    return mu_min + (mu_max - mu_min) / (1.0 + scale);
}

double PhysicsEngine::aroraMobilityHole(double T, double N) noexcept {
    const double Tr     = T / 300.0;
    const double mu_min = mu_p_min_300 * std::pow(Tr, -0.57);
    const double mu_max = mu_p_max_300 * std::pow(Tr, -2.23);
    const double N_ref  = Nref_p_300   * std::pow(Tr,  2.4);
    const double alpha  = alpha_p_300  * std::pow(Tr, -0.146);
    const double scale  = (N > 0.0) ? std::pow(N / N_ref, alpha) : 0.0;
    return mu_min + (mu_max - mu_min) / (1.0 + scale);
}


// =============================================================================
// Full-ionization solver (closed form quadratic)
// =============================================================================
void PhysicsEngine::solveFullIonization() {
    const double kT = phys::k_B * m_T;
    switch (m_dopingType) {
        case DopingType::Intrinsic:
            m_NdPlus = 0.0; m_NaMinus = 0.0;
            m_n = m_ni;     m_p = m_ni;
            break;
        case DopingType::NType: {
            const double Nd = m_N;
            m_NdPlus = Nd; m_NaMinus = 0.0;
            m_n = 0.5 * (Nd + std::sqrt(Nd * Nd + 4.0 * m_ni * m_ni));
            m_p = (m_ni * m_ni) / m_n;
            break;
        }
        case DopingType::PType: {
            const double Na = m_N;
            m_NdPlus = 0.0; m_NaMinus = Na;
            m_p = 0.5 * (Na + std::sqrt(Na * Na + 4.0 * m_ni * m_ni));
            m_n = (m_ni * m_ni) / m_p;
            break;
        }
    }
    if (m_n >= m_p) m_Ef = m_Eg - kT * std::log(m_Nc / m_n);
    else            m_Ef =        kT * std::log(m_Nv / m_p);
}


// =============================================================================
// Incomplete-ionization solver (self-consistent fixed-point on E_f)
// =============================================================================
void PhysicsEngine::solveIncompleteIonization() {
    const double kT = phys::k_B * m_T;
    const double Ec = m_Eg, Ev = 0.0;
    const double Ed = Ec - m_material->E_donor_offset;
    const double Ea = Ev + m_material->E_acceptor_offset;

    const double Nd = (m_dopingType == DopingType::NType) ? m_N : 0.0;
    const double Na = (m_dopingType == DopingType::PType) ? m_N : 0.0;

    double Ef = 0.5 * (Ec + Ev) + 0.5 * kT * std::log(m_Nv / m_Nc);

    constexpr int    maxIter = 200;
    constexpr double tol     = 1.0e-8;

    double NdPlus = Nd, NaMinus = Na, n = m_ni, p = m_ni;

    for (int iter = 0; iter < maxIter; ++iter) {
        NdPlus  = (Nd > 0.0)
                ? Nd / (1.0 + phys::g_donor    * std::exp((Ef - Ed) / kT))
                : 0.0;
        NaMinus = (Na > 0.0)
                ? Na / (1.0 + phys::g_acceptor * std::exp((Ea - Ef) / kT))
                : 0.0;

        const double dN = NdPlus - NaMinus;
        n = 0.5 * (dN + std::sqrt(dN * dN + 4.0 * m_ni * m_ni));
        p = (m_ni * m_ni) / n;

        const double Ef_new = (n >= p) ? Ec - kT * std::log(m_Nc / n)
                                       : Ev + kT * std::log(m_Nv / p);
        if (std::abs(Ef_new - Ef) < tol) { Ef = Ef_new; break; }
        Ef = 0.5 * Ef + 0.5 * Ef_new;
    }

    m_NdPlus = NdPlus; m_NaMinus = NaMinus;
    m_n = n;           m_p = p;          m_Ef = Ef;
}


// =============================================================================
// Optical generation -- direct vs indirect bandgap
//
//   Indirect (Si, Ge): alpha proportional to (hv - Eg)^2; we use a linear
//                      pre-factor here for stability over the full slider
//                      range, but with smaller K_opt.
//   Direct   (GaAs)  : alpha ~ (hv - Eg)^(1/2); much steeper onset, modelled
//                      with sqrt + larger K_opt.
//
//   The physical justification is the joint density of states:
//       g(hv) ~ (hv - Eg)^(1/2)   for direct allowed transitions
//       g(hv) ~ (hv - Eg)^2       for indirect (phonon-assisted) transitions
// =============================================================================
void PhysicsEngine::computeOptical() {
    m_Ephoton = photonEnergyEv(m_lambda_nm);

    if (!m_opticalEnabled || m_Ephoton <= m_Eg) {
        m_Gopt = m_deltaN = 0.0;
        return;
    }

    const double excess = m_Ephoton - m_Eg;
    const double power  = m_material->isDirectBandgap ? 0.5 : 1.0;
    m_deltaN = m_material->K_opt_excess * std::pow(excess, power);

    constexpr double tau_pedagogical = 1.0e-6;   // [s]
    m_Gopt = m_deltaN / tau_pedagogical;
}


// =============================================================================
// Transport: pick model, compute mobility + sigma.
//
// The Matthiessen impurity term scatters carriers off *ionized* impurities,
// so freeze-out (where N_ionized << N_d) should *increase* mobility -- and
// this fall-through here makes that happen in the simulator. The Arora model
// is left unchanged (it was empirically fitted on full-ionization data and
// has its own embedded T scaling).
// =============================================================================
void PhysicsEngine::computeTransport() {
    if (m_mobilityModel == MobilityModel::Matthiessen) {
        // Choose the relevant ionized impurity count for the current state.
        double N_scat = m_N;       // default = full ionization
        if (m_dopingType == DopingType::NType)      N_scat = m_NdPlus;
        else if (m_dopingType == DopingType::PType) N_scat = m_NaMinus;
        else                                         N_scat = 0.0;
        m_mu_n = matthiessenMobilityElectron(*m_material, m_T, N_scat);
        m_mu_p = matthiessenMobilityHole    (*m_material, m_T, N_scat);
    } else {
        m_mu_n = aroraMobilityElectron(m_T, m_N);
        m_mu_p = aroraMobilityHole    (m_T, m_N);
    }
    const double n_total = m_n + m_deltaN + m_extraDriftDN;
    const double p_total = m_p + m_deltaN + m_extraDriftDN;
    m_sigma = phys::q_e * (n_total * m_mu_n + p_total * m_mu_p);
}


// =============================================================================
// Two-carrier Hall coefficient
// =============================================================================
void PhysicsEngine::computeHall() {
    const double n = m_n + m_deltaN + m_extraDriftDN;
    const double p = m_p + m_deltaN + m_extraDriftDN;
    const double denom = phys::q_e * std::pow(p * m_mu_p + n * m_mu_n, 2);
    if (denom == 0.0) { m_R_H = 0.0; return; }
    m_R_H = (p * m_mu_p * m_mu_p - n * m_mu_n * m_mu_n) / denom;
}

double PhysicsEngine::getHallVoltage(double current_A,
                                     double thickness_cm,
                                     double field_T) const noexcept
{
    if (thickness_cm <= 0.0) return 0.0;
    return m_R_H * current_A * field_T / thickness_cm;
}


// =============================================================================
// Master update
// =============================================================================
void PhysicsEngine::recompute() {
    // Band structure (Varshni + T^(3/2) DOS scaling).
    m_Eg = bandgapAt(*m_material, m_T);
    const double T_3_2 = std::pow(m_T / 300.0, 1.5);
    m_Nc = m_material->Nc_300K * T_3_2;
    m_Nv = m_material->Nv_300K * T_3_2;
    m_ni = intrinsicCarrierAt(*m_material, m_T);

    if (m_incompleteIonization && m_dopingType != DopingType::Intrinsic)
        solveIncompleteIonization();
    else
        solveFullIonization();

    computeOptical();
    computeTransport();
    computeHall();
}


// =============================================================================
// CSV export
// =============================================================================
bool PhysicsEngine::exportCSV(const std::string& path) const {
    bool exists;
    {
        std::ifstream probe(path);
        exists = probe.good();
    }
    std::ofstream file(path, std::ios::app);
    if (!file) return false;

    if (!exists) {
        file << "Material,T[K],DopingType,N[cm^-3],IncompleteIonization,"
                "MobilityModel,Lambda[nm],OpticalOn,B[T],"
                "Eg[eV],Ef[eV],Ephoton[eV],"
                "n_i[cm^-3],n[cm^-3],p[cm^-3],Delta_n[cm^-3],DriftDN[cm^-3],"
                "Nd+[cm^-3],Na-[cm^-3],IonizationFraction,"
                "mu_n[cm^2/Vs],mu_p[cm^2/Vs],sigma[S/cm],rho[Ohm.cm],"
                "R_H[cm^3/C]\n";
    }

    const char* dopingStr =
        (m_dopingType == DopingType::Intrinsic) ? "Intrinsic" :
        (m_dopingType == DopingType::NType)     ? "NType"     : "PType";
    const char* mobilityStr =
        (m_mobilityModel == MobilityModel::Matthiessen) ? "Matthiessen"
                                                         : "Arora";

    std::ostringstream oss;
    oss << m_material->name                                       << ','
        << std::fixed      << std::setprecision(3) << m_T         << ','
        << dopingStr                                               << ','
        << std::scientific << std::setprecision(4) << m_N         << ','
        << (m_incompleteIonization ? "true" : "false")            << ','
        << mobilityStr                                             << ','
        << std::fixed      << std::setprecision(2) << m_lambda_nm << ','
        << (m_opticalEnabled ? "true" : "false")                  << ','
        << std::fixed      << std::setprecision(3) << m_B         << ','
        << std::fixed      << std::setprecision(5) << m_Eg        << ','
        << std::fixed      << std::setprecision(5) << m_Ef        << ','
        << std::fixed      << std::setprecision(5) << m_Ephoton   << ','
        << std::scientific << std::setprecision(4) << m_ni        << ','
        << std::scientific << std::setprecision(4) << m_n         << ','
        << std::scientific << std::setprecision(4) << m_p         << ','
        << std::scientific << std::setprecision(4) << m_deltaN    << ','
        << std::scientific << std::setprecision(4) << m_extraDriftDN << ','
        << std::scientific << std::setprecision(4) << m_NdPlus    << ','
        << std::scientific << std::setprecision(4) << m_NaMinus   << ','
        << std::fixed      << std::setprecision(5) << getIonizationFraction() << ','
        << std::fixed      << std::setprecision(2) << m_mu_n      << ','
        << std::fixed      << std::setprecision(2) << m_mu_p      << ','
        << std::scientific << std::setprecision(4) << m_sigma     << ','
        << std::scientific << std::setprecision(4) << getResistivity() << ','
        << std::scientific << std::setprecision(4) << m_R_H
        << '\n';

    file << oss.str();
    return file.good();
}


// =============================================================================
// Poisson / electrostatics                              [Phase 1: Poisson]
// -----------------------------------------------------------------------------
// Equilibrium Boltzmann statistics. Argument is clamped at +/-40 V_T because
// exp(40) ~ 2.4e17 is already past anything physically relevant (degenerate
// regime); IEEE-754 inf at +-700 would otherwise blow up the Gauss-Seidel.
// Reference: Sze Eq. 1.20; Pierret Eq. 4.27.
// =============================================================================
double PhysicsEngine::equilibriumElectronDensity(
    double n_i, double psi, double V_T) noexcept
{
    const double x = std::clamp(psi / V_T, -40.0, 40.0);
    return n_i * std::exp(x);
}

double PhysicsEngine::equilibriumHoleDensity(
    double n_i, double psi, double V_T) noexcept
{
    const double x = std::clamp(-psi / V_T, -40.0, 40.0);
    return n_i * std::exp(x);
}

double PhysicsEngine::builtInPotential(
    double Nd, double Na, double n_i, double V_T) noexcept
{
    if (Nd <= 0.0 || Na <= 0.0 || n_i <= 0.0) return 0.0;
    return V_T * std::log(Nd * Na / (n_i * n_i));
}

double PhysicsEngine::debyeLengthCm(
    double epsilon_r, double N_per_cm3, double V_T) noexcept
{
    if (N_per_cm3 <= 0.0) return 0.0;
    const double eps_s = phys::eps_0 * epsilon_r;
    return std::sqrt(eps_s * V_T / (phys::q_e * N_per_cm3));
}


// =============================================================================
// Non-equilibrium quasi-Fermi statistics       [Phase 3]
// -----------------------------------------------------------------------------
// Argument is clamped at +/-40 V_T (~1 V at 300 K) for numerical safety;
// this is well past degeneracy and any practical bias range.
// =============================================================================
double PhysicsEngine::nonEqElectronDensity(
    double n_i, double psi, double phi_n, double V_T) noexcept
{
    const double x = std::clamp((psi - phi_n) / V_T, -40.0, 40.0);
    return n_i * std::exp(x);
}

double PhysicsEngine::nonEqHoleDensity(
    double n_i, double psi, double phi_p, double V_T) noexcept
{
    const double x = std::clamp((phi_p - psi) / V_T, -40.0, 40.0);
    return n_i * std::exp(x);
}

// =============================================================================
// Bernoulli function (Scharfetter-Gummel core)
// -----------------------------------------------------------------------------
// Branch-clamped to keep B(x) finite over the full IEEE-754 range. For
// |x| < 1e-6 the 4-term Taylor is used (relative error < 1e-25); for
// |x| > 40 the asymptotic forms are used; in between, std::expm1 is
// preferred over (exp(x) - 1) because it preserves accuracy near zero
// even when called from optimised inlines.
// Reference: Selberherr Sec. 6.2 / Scharfetter-Gummel IEEE-ED 16 (1969).
// =============================================================================
double PhysicsEngine::bernoulli(double x) noexcept {
    constexpr double SMALL = 1.0e-6;
    constexpr double BIG   = 40.0;
    if (std::abs(x) < SMALL) return bernoulli_taylor(x);
    if (x >  BIG)            return x * std::exp(-x);  // exp(x)-1 ~ exp(x)
    if (x < -BIG)            return -x;                // exp(x) ~ 0, denom = -1
    return x / std::expm1(x);
}


// =============================================================================
// Spatially varying mobility (Phase 4)
// -----------------------------------------------------------------------------
// Matthiessen lattice + impurity scattering at local doping, then
// Caughey-Thomas saturation at the local field. Both pieces already
// exist in PhysicsEngine; this is just the canonical composition.
// =============================================================================
double PhysicsEngine::localMobilityElectron(
    const material::Profile& mat, double T,
    double N_total_per_cm3, double E_field_V_per_cm) noexcept
{
    const double mu_low = matthiessenMobilityElectron(mat, T, N_total_per_cm3);
    return highFieldMobility(mu_low, E_field_V_per_cm,
                             mat.v_sat_n, mat.beta_n);
}

double PhysicsEngine::localMobilityHole(
    const material::Profile& mat, double T,
    double N_total_per_cm3, double E_field_V_per_cm) noexcept
{
    const double mu_low = matthiessenMobilityHole(mat, T, N_total_per_cm3);
    return highFieldMobility(mu_low, E_field_V_per_cm,
                             mat.v_sat_p, mat.beta_p);
}


// =============================================================================
// Capacitance helpers (Phase 4 bonus -- small-signal C, G estimators)
// -----------------------------------------------------------------------------
// One-sided / abrupt PN junction depletion width (Sze Eq. 2.19):
//
//   W = sqrt[ 2 eps_s (V_bi - V_a) (Na + Nd) / (q Na Nd) ]
//
// Reverse bias (V_a < 0) widens W; forward bias (V_a > 0, but
// V_a < V_bi) shrinks W. We cap V_a at 0.95 V_bi to avoid the formula
// blowing up when forward bias overcomes the built-in potential
// (where the depletion approximation itself stops being valid).
// =============================================================================
double PhysicsEngine::depletionWidthFlat(
    double Nd, double Na, double V_bi, double V_a,
    double epsilon_r) noexcept
{
    if (Nd <= 0.0 || Na <= 0.0 || epsilon_r <= 0.0 || V_bi <= 0.0) return 0.0;
    const double V_eff = std::max(V_bi - V_a, 0.05 * V_bi);
    const double eps_s = phys::eps_0 * epsilon_r;
    const double Nratio = (Na + Nd) / (Na * Nd);
    return std::sqrt(2.0 * eps_s * V_eff * Nratio / phys::q_e);
}

double PhysicsEngine::depletionCapacitanceFlat(
    double Nd, double Na, double V_bi, double V_a,
    double epsilon_r, double area_cm2) noexcept
{
    const double W = depletionWidthFlat(Nd, Na, V_bi, V_a, epsilon_r);
    if (W <= 0.0) return 0.0;
    return phys::eps_0 * epsilon_r * area_cm2 / W;
}

double PhysicsEngine::diffusionCapacitance(
    double I_forward_A, double tau_s, double V_T) noexcept
{
    if (V_T <= 0.0 || I_forward_A <= 0.0 || tau_s <= 0.0) return 0.0;
    // C_d = q I tau / (kT) = I tau / V_T
    return I_forward_A * tau_s / V_T;
}


double PhysicsEngine::ohmicContactPsi(
    double V_metal, double Nd, double Na,
    double n_i, double V_T) noexcept
{
    // Charge neutrality at an ohmic contact: n - p + Na - Nd = 0,
    // combined with np = n_i^2 (thermal equilibrium at the metal). For
    // a heavily-doped contact (the practical case), this reduces to:
    //   n-contact (Nd >> Na):  psi = V_metal + V_T ln(Nd / n_i)
    //   p-contact (Na >> Nd):  psi = V_metal - V_T ln(Na / n_i)
    // Sze Eq. 2.7. We use the asymmetric form -- compensated contacts
    // are not handled (would require solving the quadratic).
    if (Nd > Na && Nd > 0.0 && n_i > 0.0) {
        return V_metal + V_T * std::log(Nd / n_i);
    }
    if (Na > 0.0 && n_i > 0.0) {
        return V_metal - V_T * std::log(Na / n_i);
    }
    return V_metal;   // intrinsic contact (uncommon): no shift.
}


// =============================================================================
// Auger recombination
// -----------------------------------------------------------------------------
// Three-particle process; the (np - n_i^2) factor enforces detailed balance
// (R = 0 in equilibrium). Dominates SRH for N > ~1e18.
// Reference: Sze Sec. 1.5.6 / Pierret Sec. 5.2.4.
// =============================================================================
double PhysicsEngine::recombAuger(
    double n, double p, double n_i,
    double C_n, double C_p) noexcept
{
    const double excess = n * p - n_i * n_i;
    return (C_n * n + C_p * p) * excess;
}


// =============================================================================
// Chynoweth impact ionization
// -----------------------------------------------------------------------------
// alpha(E) = alpha_inf * exp[-(E_crit / |E|)^m]   [cm^-1]
//
// The exponent is hard-cut at 80 (the floor of double-precision exp's useful
// range) so that vanishing E-fields produce alpha = 0 cleanly instead of
// generating subnormals or NaNs. At |E| -> 0 the original formula gives
// alpha -> 0 anyway, so the clamp does not perturb the physics.
//
// Reference: Chynoweth, Phys. Rev. 109 (1958) 1537; Sze Sec. 2.4.2;
// Van Overstraeten - de Man, Solid-State Electronics 13 (1970) 583.
// =============================================================================
double PhysicsEngine::chynowethRate(
    double E_field_V_per_cm,
    double alpha_inf, double E_crit, double m) noexcept
{
    const double E = std::abs(E_field_V_per_cm);
    if (E <= 0.0 || alpha_inf <= 0.0 || E_crit <= 0.0) return 0.0;
    const double ratio = E_crit / E;
    const double arg   = (m == 1.0) ? ratio : std::pow(ratio, m);
    if (arg > 80.0) return 0.0;        // exp underflow guard
    return alpha_inf * std::exp(-arg);
}


// =============================================================================
// Avalanche multiplication factor (uniform field, single-carrier limit)
// -----------------------------------------------------------------------------
// For arbitrary alpha_n != alpha_p the exact ionization integral gives
//
//     1 - 1/M = int_0^W alpha_n exp{ -int_0^x (alpha_n - alpha_p) dx' } dx
//
// In the homogeneous-field, single-carrier limit (alpha_n ~ alpha_p =
// alpha_eff) this collapses to the classical Miller form
//
//     M = 1 / (1 - alpha_eff * W)
//
// which we use here because it is the tractable closed-form expression
// suitable for a real-time UI gauge.  alpha_eff is taken as max(alpha_n,
// alpha_p) so the answer is conservative (breakdown predicted earlier).
//
// Reference: Sze Sec. 2.4.2 Eq. 65; Miller, Phys. Rev. 99 (1955) 1234.
// =============================================================================
double PhysicsEngine::avalancheMultiplication(
    double alpha_n, double alpha_p,
    double depletion_width_cm) noexcept
{
    const double a_eff = std::max(alpha_n, alpha_p);
    const double prod  = a_eff * std::max(0.0, depletion_width_cm);
    if (prod >= 0.999) return std::numeric_limits<double>::infinity();
    return 1.0 / (1.0 - prod);
}


// =============================================================================
// Anderson rule -- valence-band edge (Phase 5 helper)
// -----------------------------------------------------------------------------
//   E_v(x) = E_c(x) - E_g(x)
// where E_c(x) is the static helper above (in absolute eV with the engine's
// reference material as the chi anchor). E_g is evaluated at the requested
// temperature so hot cells display narrowed gaps in BandView. The chi anchor
// is the reference profile, so the *relative* step at A/B interfaces equals
//   Delta E_v(B - A) = (chi_A - chi_B) + (E_g,A - E_g,B)
// = -(Delta E_c) - (Delta E_g of A->B) -- matching andersonOffsets().
// =============================================================================
double PhysicsEngine::valenceBandEdge(
    double psi_volts,
    const material::Profile& mat_local,
    const material::Profile& mat_ref,
    double T) noexcept
{
    const double Ec = conductionBandEdge(psi_volts, mat_local, mat_ref);
    const double Eg = bandgapAt(mat_local, T);
    return Ec - Eg;
}


// =============================================================================
// Band-to-band (Zener) tunneling -- Kane model
// -----------------------------------------------------------------------------
//
//   G_BTBT(E) = A * E^P * exp(-B / |E|)         [cm^-3 / s]
//
//     P = 2     (direct gap)
//     P = 5/2   (indirect, phonon-assisted)
//
// At 300 K in Si this becomes appreciable above |E| ~ 1e6 V/cm; in Ge or
// narrow-gap III-V it lights up an order of magnitude earlier (smaller B).
//
// Reference: Kane, J. Phys. Chem. Solids 12 (1959) 181;
//            Hurkx et al., IEEE Trans. ED 39 (1992) 331; Sze Sec. 4.3.
// =============================================================================
double PhysicsEngine::kaneBTBT(
    double E_field_V_per_cm,
    double A_kane, double B_kane,
    bool isDirect) noexcept
{
    const double E = std::abs(E_field_V_per_cm);
    if (E <= 0.0 || A_kane <= 0.0 || B_kane <= 0.0) return 0.0;
    const double exp_arg = B_kane / E;
    if (exp_arg > 80.0) return 0.0;       // tunneling negligible
    const double power = isDirect ? (E * E)
                                  : (E * E * std::sqrt(E));   // E^(5/2)
    return A_kane * power * std::exp(-exp_arg);
}
