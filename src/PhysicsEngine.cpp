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
// Transport: pick model, compute mobility + sigma
// =============================================================================
void PhysicsEngine::computeTransport() {
    if (m_mobilityModel == MobilityModel::Matthiessen) {
        m_mu_n = matthiessenMobilityElectron(*m_material, m_T, m_N);
        m_mu_p = matthiessenMobilityHole    (*m_material, m_T, m_N);
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
    // Bandgap (Varshni, parameterized by material).
    m_Eg = m_material->Eg0
         - (m_material->varshni_a * m_T * m_T) / (m_T + m_material->varshni_b);

    const double T_3_2 = std::pow(m_T / 300.0, 1.5);
    m_Nc = m_material->Nc_300K * T_3_2;
    m_Nv = m_material->Nv_300K * T_3_2;

    const double kT = phys::k_B * m_T;
    m_ni = std::sqrt(m_Nc * m_Nv) * std::exp(-m_Eg / (2.0 * kT));

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
