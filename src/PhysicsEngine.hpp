#pragma once

// =============================================================================
// PhysicsEngine.hpp -- multi-material semiconductor physics core
//
//   Author : dex / cemfly-april2026
//   License: MIT
// -----------------------------------------------------------------------------
// Semiconductor physics core for the Analysis & Simulation Platform.
//
// All material-dependent parameters are now sourced from a Material::Profile
// (see Material.hpp).  The engine itself only encodes the *physics* and
// statistics; it can be retargeted to Si, GaAs or Ge at runtime via
// setMaterial().
//
// References
//   [1] Kasap, "Principles of Electronic Materials and Devices", 4th ed.
//   [2] Kittel, "Introduction to Solid State Physics", 8th ed.
//   [3] Sze & Ng, "Physics of Semiconductor Devices", 3rd ed.
// =============================================================================

#include <cmath>
#include <string>

#include "Material.hpp"


// =============================================================================
// Universal physical constants and UI clamp ranges (material-independent).
// =============================================================================
namespace phys {

    inline constexpr double k_B     = 8.617333262e-5;   // [eV/K]
    inline constexpr double q_e     = 1.602176634e-19;  // [C]
    inline constexpr double h_eVs   = 4.135667696e-15;  // [eV s]
    inline constexpr double c_light = 2.99792458e8;     // [m/s]
    inline constexpr double hc_eVnm = 1239.84198;       // E[eV] = hc/lambda[nm]

    inline constexpr double g_donor    = 2.0;           // donor degeneracy
    inline constexpr double g_acceptor = 4.0;           // acceptor degeneracy

    inline constexpr double lambda_min_nm =  300.0;
    inline constexpr double lambda_max_nm = 1500.0;
    inline constexpr double B_min_T       = 0.0;
    inline constexpr double B_max_T       = 10.0;

} // namespace phys


// =============================================================================
// Doping regime selected by the UI.
// =============================================================================
enum class DopingType {
    Intrinsic,
    NType,
    PType
};


// =============================================================================
// Mobility model selector.
// =============================================================================
enum class MobilityModel {
    Matthiessen,   // 1/mu = 1/mu_L(T) + 1/mu_I(T,N)  -- physics-grounded
    Arora          // empirical Caughey-Thomas/Arora (Si-tuned, legacy)
};


// =============================================================================
// PhysicsEngine
// =============================================================================
class PhysicsEngine {
public:
    explicit PhysicsEngine(material::Kind kind = material::Kind::Silicon);

    // ---- User inputs -------------------------------------------------------
    void setMaterial             (material::Kind kind);
    void setTemperature          (double T_kelvin);
    void setDopingType           (DopingType type);
    void setDopingConcentration  (double N_per_cm3);
    void setIncompleteIonization (bool enabled);
    void setMobilityModel        (MobilityModel m);
    void setWavelengthNm         (double lambda_nm);
    void setOpticalEnabled       (bool enabled);
    void setMagneticField        (double B_tesla);

    // Drift-diffusion: an externally-managed source bumps the steady-state
    // excess carrier density. Updated by the DriftDiffusion module each frame.
    void setDriftDiffusionExcess (double extra_dN);

    // ---- Query: user inputs ------------------------------------------------
    [[nodiscard]] material::Kind getMaterialKind()      const noexcept { return m_materialKind;       }
    [[nodiscard]] const material::Profile& getMaterial()const noexcept { return *m_material;          }
    [[nodiscard]] double         getTemperature()       const noexcept { return m_T;                  }
    [[nodiscard]] DopingType     getDopingType()        const noexcept { return m_dopingType;         }
    [[nodiscard]] double         getDopingConcentration()const noexcept{ return m_N;                  }
    [[nodiscard]] bool           getIncompleteIonization()const noexcept{return m_incompleteIonization;}
    [[nodiscard]] MobilityModel  getMobilityModel()     const noexcept { return m_mobilityModel;      }
    [[nodiscard]] double         getWavelengthNm()      const noexcept { return m_lambda_nm;          }
    [[nodiscard]] bool           getOpticalEnabled()    const noexcept { return m_opticalEnabled;     }
    [[nodiscard]] double         getMagneticField()     const noexcept { return m_B;                  }
    [[nodiscard]] double         getDriftDiffusionExcess()const noexcept{return m_extraDriftDN;       }

    // ---- Query: thermodynamic state ----------------------------------------
    [[nodiscard]] double getBandgap()               const noexcept { return m_Eg; }
    [[nodiscard]] double getEffectiveNc()           const noexcept { return m_Nc; }
    [[nodiscard]] double getEffectiveNv()           const noexcept { return m_Nv; }
    [[nodiscard]] double getIntrinsicCarrier()      const noexcept { return m_ni; }
    [[nodiscard]] double getElectronConcentration() const noexcept { return m_n;  }
    [[nodiscard]] double getHoleConcentration()     const noexcept { return m_p;  }
    [[nodiscard]] double getFermiLevel()            const noexcept { return m_Ef; }

    // ---- Query: ionization -------------------------------------------------
    [[nodiscard]] double getIonizedDonors()      const noexcept { return m_NdPlus;  }
    [[nodiscard]] double getIonizedAcceptors()   const noexcept { return m_NaMinus; }
    [[nodiscard]] double getIonizationFraction() const noexcept;

    // ---- Query: transport --------------------------------------------------
    [[nodiscard]] double getElectronMobility() const noexcept { return m_mu_n;  }
    [[nodiscard]] double getHoleMobility()     const noexcept { return m_mu_p;  }
    [[nodiscard]] double getConductivity()     const noexcept { return m_sigma; }
    [[nodiscard]] double getResistivity()      const noexcept;

    [[nodiscard]] double getTotalElectronConc() const noexcept { return m_n + m_deltaN + m_extraDriftDN; }
    [[nodiscard]] double getTotalHoleConc()     const noexcept { return m_p + m_deltaN + m_extraDriftDN; }

    // ---- Query: optical ----------------------------------------------------
    [[nodiscard]] double getPhotonEnergy()         const noexcept { return m_Ephoton; }
    [[nodiscard]] double getOpticalGeneration()    const noexcept { return m_Gopt;    }
    [[nodiscard]] double getExcessCarrierDensity() const noexcept { return m_deltaN;  }
    [[nodiscard]] bool   isOpticallyPumped()       const noexcept;

    // ---- Query: Hall / magnetic --------------------------------------------
    [[nodiscard]] double getHallCoefficient() const noexcept { return m_R_H; }
    [[nodiscard]] double getHallVoltage(double current_A,
                                        double thickness_cm,
                                        double field_T) const noexcept;

    // ---- Band / dopant levels ---------------------------------------------
    [[nodiscard]] double getValenceBandEdge()    const noexcept { return 0.0;  }
    [[nodiscard]] double getConductionBandEdge() const noexcept { return m_Eg; }
    [[nodiscard]] double getDonorLevel()         const noexcept;
    [[nodiscard]] double getAcceptorLevel()      const noexcept;

    [[nodiscard]] double fermiDirac(double E) const noexcept;

    // ---- CSV export --------------------------------------------------------
    [[nodiscard]] bool exportCSV(const std::string& path) const;

    // ---- Static helpers (parameterized by material; testable directly) ----
    //
    // Numerically robust evaluators -- see PhysicsEngine.cpp for the
    // log-domain implementation that survives T < 50 K (where the simple
    // exp(-Eg/2kT) form underflows to zero in IEEE754).
    [[nodiscard]] static double bandgapAt(
        const material::Profile& mat, double T) noexcept;
    [[nodiscard]] static double intrinsicCarrierAt(
        const material::Profile& mat, double T) noexcept;
    [[nodiscard]] static double log10IntrinsicCarrierAt(
        const material::Profile& mat, double T) noexcept;

    // The N argument is the *ionized* impurity concentration -- which is
    // what actually scatters carriers. At freeze-out (low T) this is much
    // smaller than the dopant concentration, and the resulting mobility is
    // correspondingly higher.
    [[nodiscard]] static double matthiessenMobilityElectron(
        const material::Profile& mat, double T, double N_ionized) noexcept;
    [[nodiscard]] static double matthiessenMobilityHole(
        const material::Profile& mat, double T, double N_ionized) noexcept;
    [[nodiscard]] static double aroraMobilityElectron(double T, double N) noexcept;
    [[nodiscard]] static double aroraMobilityHole    (double T, double N) noexcept;

    // Caughey-Thomas high-field saturation:
    //   v(E) = mu_low * E / [1 + (mu_low * E / v_sat)^beta]^(1/beta)
    // returned as the effective field-dependent mobility v(E)/E.
    [[nodiscard]] static double highFieldMobility(
        double mu_low, double E_field_V_per_cm,
        double v_sat_cm_per_s, double beta) noexcept;

    // Shockley-Read-Hall recombination through midgap traps:
    //   R = (n p - n_i^2) / [tau_p (n + n_i) + tau_n (p + n_i)]
    // (units: [cm^-3 / s])
    [[nodiscard]] static double recombSRH(
        double n, double p, double n_i,
        double tau_n, double tau_p) noexcept;

    // Slotboom-de Graaff bandgap narrowing -- significant for N > ~1e18.
    //   dE_g = 9e-3 * [ln(N/N_ref) + sqrt(ln(N/N_ref)^2 + 0.5)]   eV
    [[nodiscard]] static double bandgapNarrowing(double N_per_cm3) noexcept;

    // Effective intrinsic carrier under heavy doping:  n_ie = n_i * exp(dE_g/2kT).
    [[nodiscard]] static double effectiveIntrinsicCarrier(
        const material::Profile& mat, double T, double N_per_cm3) noexcept;

    // ---- Optical: Tauc-style absorption coefficient and penetration depth -
    //
    //   alpha(hv) = A * (hv - E_g)^p / hv         (p = 2 indirect, 1/2 direct)
    //   penetration depth L_alpha = 1 / alpha
    //
    //   Reference: Pankove, "Optical Processes in Semiconductors" Ch. 3;
    //   Sze App. C.
    [[nodiscard]] static double absorptionCoefficient(
        const material::Profile& mat, double photon_eV) noexcept;
    [[nodiscard]] static double penetrationDepthCm(
        const material::Profile& mat, double photon_eV) noexcept;

    // ---- BJT (NPN, textbook one-sided emitter junction) -------------------
    //
    //   gamma  = 1 / (1 + (D_p N_B W_E) / (D_n N_E W_B))     (Sze 5.2)
    //   alpha_T = 1 - W_B^2 / (2 L_n^2)                       (transport factor)
    //   beta   = gamma * alpha_T / (1 - gamma * alpha_T)
    //   I_C(V_CE) = I_C0 * (1 + V_CE / V_A)                   (Early effect)
    [[nodiscard]] static double bjtEmitterEfficiency(
        double D_n, double N_E, double W_E,
        double D_p, double N_B, double W_B) noexcept;
    [[nodiscard]] static double bjtCurrentGain(
        double gamma, double alpha_T) noexcept;
    [[nodiscard]] static double earlyEffectFactor(
        double V_CE, double V_Early) noexcept;

    [[nodiscard]] static double photonEnergyEv(double lambda_nm) noexcept;

private:
    void recompute();
    void solveFullIonization();
    void solveIncompleteIonization();
    void computeOptical();
    void computeTransport();
    void computeHall();

    // ---- Inputs ------------------------------------------------------------
    material::Kind            m_materialKind;
    const material::Profile*  m_material;

    double        m_T;
    DopingType    m_dopingType;
    double        m_N;
    bool          m_incompleteIonization;
    MobilityModel m_mobilityModel;

    double        m_lambda_nm;
    bool          m_opticalEnabled;
    double        m_B;
    double        m_extraDriftDN = 0.0;     // injected by DriftDiffusion

    // ---- Outputs -----------------------------------------------------------
    double m_Eg, m_Nc, m_Nv, m_ni;
    double m_n,  m_p,  m_Ef;
    double m_NdPlus, m_NaMinus;
    double m_mu_n, m_mu_p, m_sigma;
    double m_Ephoton, m_Gopt, m_deltaN;
    double m_R_H;
};
