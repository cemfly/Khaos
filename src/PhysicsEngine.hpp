#pragma once

// =============================================================================
// PhysicsEngine.hpp -- Phase 6: multi-material core
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

    // ---- Static mobility helpers (parameterized by material) --------------
    [[nodiscard]] static double matthiessenMobilityElectron(
        const material::Profile& mat, double T, double N) noexcept;
    [[nodiscard]] static double matthiessenMobilityHole(
        const material::Profile& mat, double T, double N) noexcept;
    [[nodiscard]] static double aroraMobilityElectron(double T, double N) noexcept;
    [[nodiscard]] static double aroraMobilityHole    (double T, double N) noexcept;

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
