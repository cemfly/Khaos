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
PhysicsEngine::PhysicsEngine()
    : m_T(300.0)
    , m_dopingType(DopingType::Intrinsic)
    , m_N(0.0)
    , m_incompleteIonization(false)
    , m_mobilityModel(MobilityModel::Matthiessen)
    , m_lambda_nm(1500.0)          // start below bandgap => no illumination
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
// Setters: every input change triggers a full pipeline recompute so the
// cached outputs stay internally consistent.
// =============================================================================
void PhysicsEngine::setTemperature(double T) {
    m_T = std::max(1.0, T);          // keep kT finite
    recompute();
}

void PhysicsEngine::setDopingType(DopingType t) {
    m_dopingType = t;
    recompute();
}

void PhysicsEngine::setDopingConcentration(double N) {
    m_N = std::max(0.0, N);
    recompute();
}

void PhysicsEngine::setIncompleteIonization(bool e) {
    m_incompleteIonization = e;
    recompute();
}

void PhysicsEngine::setMobilityModel(MobilityModel m) {
    m_mobilityModel = m;
    recompute();
}

void PhysicsEngine::setWavelengthNm(double lambda) {
    m_lambda_nm = std::clamp(lambda, phys::lambda_min_nm, phys::lambda_max_nm);
    recompute();
}

void PhysicsEngine::setOpticalEnabled(bool e) {
    m_opticalEnabled = e;
    recompute();
}

void PhysicsEngine::setMagneticField(double B) {
    m_B = std::clamp(B, phys::B_min_T, phys::B_max_T);
    recompute();
}


// =============================================================================
// Fermi-Dirac distribution
//     f(E) = 1 / (1 + exp((E - E_f) / kT))                      [Kittel Ch. 6]
// =============================================================================
double PhysicsEngine::fermiDirac(double E) const {
    const double kT = phys::k_B * m_T;
    const double x  = (E - m_Ef) / kT;
    if (x >  50.0) return 0.0;
    if (x < -50.0) return 1.0;
    return 1.0 / (1.0 + std::exp(x));
}


// =============================================================================
// Ionization fraction  (N_d+/N_d for n-type, N_a-/N_a for p-type).
// =============================================================================
double PhysicsEngine::getIonizationFraction() const {
    if (m_dopingType == DopingType::Intrinsic || m_N <= 0.0) return 1.0;
    if (m_dopingType == DopingType::NType) return m_NdPlus  / m_N;
    return                                        m_NaMinus / m_N;
}

double PhysicsEngine::getResistivity() const {
    return (m_sigma > 0.0) ? 1.0 / m_sigma : std::numeric_limits<double>::infinity();
}

bool PhysicsEngine::isOpticallyPumped() const {
    return m_opticalEnabled && (m_Ephoton > m_Eg);
}


// =============================================================================
// Static helper: photon energy from wavelength
//     E[eV] = h c / lambda    -> with lambda in nm:  E = 1239.84 / lambda
// =============================================================================
double PhysicsEngine::photonEnergyEv(double lambda_nm) {
    return (lambda_nm > 0.0) ? phys::hc_eVnm / lambda_nm : 0.0;
}


// =============================================================================
// Matthiessen mobility model.
//
//   The key physical idea (Kittel Ch. 6.9): independent scattering mechanisms
//   add as inverse mobilities,
//
//       1 / mu_total = 1 / mu_lattice(T) + 1 / mu_impurity(T, N)
//
//   Each mechanism has its own temperature dependence:
//
//       mu_L(T)   propto T^(-3/2)    -- acoustic-phonon scattering
//                                       (Kasap 2.6, Sze Ch. 1)
//       mu_I(T,N) propto T^( 3/2)/N  -- Conwell-Weisskopf ionized-impurity
//                                       scattering  (Kittel, Sze)
//
//   At low T and low N the lattice term dominates (mu large, grows as T drops
//   until impurities freeze out). At low T and high N the impurity term
//   becomes limiting -- this is why doped silicon does not become arbitrarily
//   conductive at cryogenic temperatures.
// =============================================================================
double PhysicsEngine::matthiessenMobilityElectron(double T, double N) {
    const double Tr   = T / 300.0;
    const double mu_L = phys::mu_L_n_300 * std::pow(Tr, -1.5);

    if (N <= 0.0) return mu_L;          // intrinsic: impurity scattering absent

    const double mu_I = phys::mu_I_n_300
                      * std::pow(Tr, 1.5)
                      * (phys::N_ref_matt / N);
    return 1.0 / (1.0 / mu_L + 1.0 / mu_I);
}

double PhysicsEngine::matthiessenMobilityHole(double T, double N) {
    const double Tr   = T / 300.0;
    const double mu_L = phys::mu_L_p_300 * std::pow(Tr, -1.5);

    if (N <= 0.0) return mu_L;

    const double mu_I = phys::mu_I_p_300
                      * std::pow(Tr, 1.5)
                      * (phys::N_ref_matt / N);
    return 1.0 / (1.0 / mu_L + 1.0 / mu_I);
}


// =============================================================================
// Arora empirical model (kept for comparison with Matthiessen):
//   mu(N,T) = mu_min(T) + (mu_max(T) - mu_min(T)) / (1 + (N/N_ref(T))^alpha)
//   with T-scaling exponents from the Arora-Hauser-Roulston 1982 paper.
// =============================================================================
double PhysicsEngine::aroraMobilityElectron(double T, double N) {
    const double Tr     = T / 300.0;
    const double mu_min = phys::mu_n_min_300 * std::pow(Tr, -0.57);
    const double mu_max = phys::mu_n_max_300 * std::pow(Tr, -2.33);
    const double N_ref  = phys::Nref_n_300   * std::pow(Tr,  2.4);
    const double alpha  = phys::alpha_n_300  * std::pow(Tr, -0.146);
    const double scale  = (N > 0.0) ? std::pow(N / N_ref, alpha) : 0.0;
    return mu_min + (mu_max - mu_min) / (1.0 + scale);
}

double PhysicsEngine::aroraMobilityHole(double T, double N) {
    const double Tr     = T / 300.0;
    const double mu_min = phys::mu_p_min_300 * std::pow(Tr, -0.57);
    const double mu_max = phys::mu_p_max_300 * std::pow(Tr, -2.23);
    const double N_ref  = phys::Nref_p_300   * std::pow(Tr,  2.4);
    const double alpha  = phys::alpha_p_300  * std::pow(Tr, -0.146);
    const double scale  = (N > 0.0) ? std::pow(N / N_ref, alpha) : 0.0;
    return mu_min + (mu_max - mu_min) / (1.0 + scale);
}


// =============================================================================
// Full-ionization solver (closed form).
//   Charge neutrality  n - p = N_d - N_a   combined with mass action law
//   n p = n_i^2   yields a single quadratic for the majority carrier.
// =============================================================================
void PhysicsEngine::solveFullIonization() {
    const double kT = phys::k_B * m_T;

    switch (m_dopingType) {
        case DopingType::Intrinsic:
            m_NdPlus  = 0.0;  m_NaMinus = 0.0;
            m_n = m_ni;       m_p = m_ni;
            break;
        case DopingType::NType: {
            const double Nd = m_N;
            m_NdPlus  = Nd;   m_NaMinus = 0.0;
            m_n = 0.5 * (Nd + std::sqrt(Nd * Nd + 4.0 * m_ni * m_ni));
            m_p = (m_ni * m_ni) / m_n;
            break;
        }
        case DopingType::PType: {
            const double Na = m_N;
            m_NdPlus  = 0.0;  m_NaMinus = Na;
            m_p = 0.5 * (Na + std::sqrt(Na * Na + 4.0 * m_ni * m_ni));
            m_n = (m_ni * m_ni) / m_p;
            break;
        }
    }

    // Majority-branch Fermi level (better conditioned than minority branch).
    if (m_n >= m_p) m_Ef = m_Eg - kT * std::log(m_Nc / m_n);
    else            m_Ef =        kT * std::log(m_Nv / m_p);
}


// =============================================================================
// Incomplete-ionization solver (self-consistent).
//
//   Donor ionization:    N_d^+ = N_d / (1 + g_D exp((E_f - E_d)/kT))
//   Acceptor ionization: N_a^- = N_a / (1 + g_A exp((E_a - E_f)/kT))
//                                                                  [Kasap 5.8]
//
//   Charge neutrality:   n + N_a^- = p + N_d^+
//   Mass action law:     n p = n_i^2
//
//   Because N_d^+, N_a^- depend on E_f, the system is transcendental and
//   must be solved iteratively.  We use a damped fixed-point iteration
//   seeded with the intrinsic Fermi level; convergence is reached in O(10)
//   iterations across the whole T-N parameter space.
// =============================================================================
void PhysicsEngine::solveIncompleteIonization() {
    const double kT = phys::k_B * m_T;
    const double Ec = m_Eg;
    const double Ev = 0.0;
    const double Ed = Ec - phys::E_donor_offset;
    const double Ea = Ev + phys::E_acceptor_offset;

    const double Nd = (m_dopingType == DopingType::NType) ? m_N : 0.0;
    const double Na = (m_dopingType == DopingType::PType) ? m_N : 0.0;

    double Ef = 0.5 * (Ec + Ev) + 0.5 * kT * std::log(m_Nv / m_Nc);

    constexpr int    maxIter = 200;
    constexpr double tol     = 1.0e-8;

    double NdPlus  = Nd;
    double NaMinus = Na;
    double n = m_ni, p = m_ni;

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
        Ef = 0.5 * Ef + 0.5 * Ef_new;       // damped update
    }

    m_NdPlus  = NdPlus;  m_NaMinus = NaMinus;
    m_n       = n;       m_p       = p;     m_Ef = Ef;
}


// =============================================================================
// Optical generation / excess-carrier model.
//
//   E_photon = hc / lambda                                     [Kittel Ch. 15]
//
//   Interband absorption requires hv > E_g.  Pedagogically, we bundle the
//   photon flux, absorption coefficient and recombination lifetime into a
//   single scaling factor K_opt so that the steady-state excess is
//
//       Delta n = Delta p = K_opt * (hv - E_g)   for  hv > E_g
//
//   Both electrons and holes are generated in pairs, so the injected
//   concentration is the same.
// =============================================================================
void PhysicsEngine::computeOptical() {
    m_Ephoton = photonEnergyEv(m_lambda_nm);

    if (!m_opticalEnabled || m_Ephoton <= m_Eg) {
        m_Gopt   = 0.0;
        m_deltaN = 0.0;
        return;
    }

    const double excess_eV = m_Ephoton - m_Eg;
    m_deltaN = phys::K_opt_excess * excess_eV;

    // Purely informational: "pseudo generation rate" = Delta n / tau_arbitrary.
    // We pick tau = 1 us so G_opt reads out in units physicists recognise.
    constexpr double tau_pedagogical = 1.0e-6;   // [s]
    m_Gopt = m_deltaN / tau_pedagogical;          // [cm^-3/s]
}


// =============================================================================
// Transport: mobility dispatch, conductivity, resistivity.
//
//   sigma = q (n_total mu_n + p_total mu_p)      n_total = n + Delta n, etc.
//
//   Under illumination the excess carriers bump sigma up noticeably -- this
//   is the physical basis of a photoconductor or a solar-cell short-circuit.
// =============================================================================
void PhysicsEngine::computeTransport() {
    if (m_mobilityModel == MobilityModel::Matthiessen) {
        m_mu_n = matthiessenMobilityElectron(m_T, m_N);
        m_mu_p = matthiessenMobilityHole    (m_T, m_N);
    } else {
        m_mu_n = aroraMobilityElectron(m_T, m_N);
        m_mu_p = aroraMobilityHole    (m_T, m_N);
    }

    const double n_total = m_n + m_deltaN;
    const double p_total = m_p + m_deltaN;
    m_sigma = phys::q_e * (n_total * m_mu_n + p_total * m_mu_p);
}


// =============================================================================
// Hall coefficient (two-carrier formula).
//
//   R_H = (p mu_p^2 - n mu_n^2) / (q (p mu_p + n mu_n)^2)       [Kasap 2.4]
//
//   Sign convention:
//       R_H < 0    => majority electrons (n-type)
//       R_H > 0    => majority holes     (p-type)
//
//   Units:  cm^3 / C  (SI would be m^3/C; we keep cgs-cm for consistency
//           with the rest of the simulator which uses cm for concentrations).
// =============================================================================
void PhysicsEngine::computeHall() {
    const double n = m_n + m_deltaN;
    const double p = m_p + m_deltaN;
    const double denom = phys::q_e * std::pow(p * m_mu_p + n * m_mu_n, 2);
    if (denom == 0.0) { m_R_H = 0.0; return; }
    m_R_H = (p * m_mu_p * m_mu_p - n * m_mu_n * m_mu_n) / denom;
}

double PhysicsEngine::getHallVoltage(double I_A,
                                     double thickness_cm,
                                     double B_T) const
{
    // V_H = R_H * I * B / t       (Kasap Eq. 2.49)
    if (thickness_cm <= 0.0) return 0.0;
    return m_R_H * I_A * B_T / thickness_cm;
}


// =============================================================================
// Master update.
// =============================================================================
void PhysicsEngine::recompute() {
    // ---- Temperature-dependent band structure -----------------------------
    m_Eg = phys::Eg0_Si
         - (phys::varshni_a * m_T * m_T) / (m_T + phys::varshni_b);

    const double T_3_2 = std::pow(m_T / 300.0, 1.5);
    m_Nc = phys::Nc_300K * T_3_2;
    m_Nv = phys::Nv_300K * T_3_2;

    const double kT = phys::k_B * m_T;
    m_ni = std::sqrt(m_Nc * m_Nv) * std::exp(-m_Eg / (2.0 * kT));

    // ---- Carriers + Fermi level -------------------------------------------
    if (m_incompleteIonization && m_dopingType != DopingType::Intrinsic) {
        solveIncompleteIonization();
    } else {
        solveFullIonization();
    }

    // ---- Optical generation (uses m_Eg) -----------------------------------
    computeOptical();

    // ---- Transport (uses m_n, m_p, m_deltaN) ------------------------------
    computeTransport();

    // ---- Hall coefficient (uses transport output) -------------------------
    computeHall();
}


// =============================================================================
// CSV export.
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
        file << "T[K],"
                "DopingType,"
                "N[cm^-3],"
                "IncompleteIonization,"
                "MobilityModel,"
                "Lambda[nm],"
                "OpticalOn,"
                "B[T],"
                "Eg[eV],"
                "Ef[eV],"
                "Ephoton[eV],"
                "n_i[cm^-3],"
                "n[cm^-3],"
                "p[cm^-3],"
                "Delta_n[cm^-3],"
                "Nd+[cm^-3],"
                "Na-[cm^-3],"
                "IonizationFraction,"
                "mu_n[cm^2/Vs],"
                "mu_p[cm^2/Vs],"
                "sigma[S/cm],"
                "rho[Ohm.cm],"
                "R_H[cm^3/C]\n";
    }

    const char* dopingStr =
        (m_dopingType == DopingType::Intrinsic) ? "Intrinsic" :
        (m_dopingType == DopingType::NType)     ? "NType"     : "PType";

    const char* mobilityStr =
        (m_mobilityModel == MobilityModel::Matthiessen) ? "Matthiessen"
                                                         : "Arora";

    std::ostringstream oss;
    oss << std::fixed      << std::setprecision(3) << m_T           << ','
        << dopingStr                                                 << ','
        << std::scientific << std::setprecision(4) << m_N           << ','
        << (m_incompleteIonization ? "true" : "false")              << ','
        << mobilityStr                                               << ','
        << std::fixed      << std::setprecision(2) << m_lambda_nm   << ','
        << (m_opticalEnabled ? "true" : "false")                    << ','
        << std::fixed      << std::setprecision(3) << m_B           << ','
        << std::fixed      << std::setprecision(5) << m_Eg          << ','
        << std::fixed      << std::setprecision(5) << m_Ef          << ','
        << std::fixed      << std::setprecision(5) << m_Ephoton     << ','
        << std::scientific << std::setprecision(4) << m_ni          << ','
        << std::scientific << std::setprecision(4) << m_n           << ','
        << std::scientific << std::setprecision(4) << m_p           << ','
        << std::scientific << std::setprecision(4) << m_deltaN      << ','
        << std::scientific << std::setprecision(4) << m_NdPlus      << ','
        << std::scientific << std::setprecision(4) << m_NaMinus     << ','
        << std::fixed      << std::setprecision(5) << getIonizationFraction() << ','
        << std::fixed      << std::setprecision(2) << m_mu_n        << ','
        << std::fixed      << std::setprecision(2) << m_mu_p        << ','
        << std::scientific << std::setprecision(4) << m_sigma       << ','
        << std::scientific << std::setprecision(4) << getResistivity() << ','
        << std::scientific << std::setprecision(4) << m_R_H
        << '\n';

    file << oss.str();
    return file.good();
}
