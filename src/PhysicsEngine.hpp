#pragma once

// =============================================================================
// PhysicsEngine.hpp    (Phase 3: multidisciplinary solid-state simulator)
// -----------------------------------------------------------------------------
// Semiconductor physics core for the Interactive Bandgap & Doping Simulator.
//
// References
// ----------
//   [1] S. O. Kasap, "Principles of Electronic Materials and Devices",
//       4th ed., McGraw-Hill, 2017.
//   [2] C. Kittel, "Introduction to Solid State Physics", 8th ed.,
//       Wiley, 2005.
//   [3] Arora, Hauser, Roulston, IEEE TED-29 (1982) -- Si mobility.
//   [4] Sze & Ng, "Physics of Semiconductor Devices", 3rd ed., Wiley.
//
// Phase 1   (thermal)
// -------------------
//   * Varshni bandgap E_g(T)                                    [Kasap 5.5]
//   * Effective DOS  N_c(T), N_v(T)  with T^(3/2) scaling        [Kasap 5.1]
//   * Intrinsic carrier n_i(T)                                   [Kasap 5.2]
//   * Charge neutrality + Mass Action Law                        [Kasap 5.6]
//   * Fermi-Dirac distribution                                   [Kittel 6]
//
// Phase 2   (extended transport + reliability)
// --------------------------------------------
//   * Incomplete ionization (freeze-out)                         [Kasap 5.8]
//   * Arora empirical mobility  mu(N, T)                         [Ref. 3]
//   * Conductivity  sigma = q (n mu_n + p mu_p)                  [Kasap 2.2]
//   * CSV data export of every thermodynamic / electrical state
//
// Phase 3   (multidisciplinary: optical, magnetic, transport refinement)
// ----------------------------------------------------------------------
//   * Optical absorption:
//         E_photon = h c / lambda                                 [Kittel 15]
//         if E_photon > E_g  => electron-hole pair generation
//         Steady-state excess  Delta n = G_opt * tau
//   * Magnetic field & Hall effect:
//         F_Lorentz = q (v x B)                                   [Kasap 2.4]
//         Hall coefficient R_H = (p mu_p^2 - n mu_n^2) / (q (p mu_p + n mu_n)^2)
//   * Matthiessen's rule for mobility (explicit decomposition):
//         1 / mu = 1 / mu_lattice + 1 / mu_impurity              [Kittel 6.9]
//       with
//         mu_L        ~ T^(-3/2)       (acoustic-phonon scattering)
//         mu_I        ~ T^( 3/2) / N   (Conwell-Weisskopf impurity
//                                       scattering)
//     Arora model kept as a static helper for comparison / tests.
//
// Energy reference: E_v = 0, E_c = E_g.
// =============================================================================

#include <cmath>
#include <string>


// =============================================================================
// Physical constants & material parameters.
// =============================================================================
namespace phys {

    // ---- Fundamental constants ---------------------------------------------
    inline constexpr double k_B     = 8.617333262e-5;   // Boltzmann       [eV/K]
    inline constexpr double q_e     = 1.602176634e-19;  // Elementary charge [C]
    inline constexpr double h_eVs   = 4.135667696e-15;  // Planck h         [eV s]
    inline constexpr double c_light = 2.99792458e8;     // Speed of light   [m/s]

    // Handy convenience: hc expressed in eV.nm so  E_photon[eV] = hc_eVnm/lambda[nm]
    inline constexpr double hc_eVnm = 1239.84198;

    // ---- Silicon bulk parameters (room-temperature reference) --------------
    inline constexpr double Nc_300K = 2.80e19;   // [cm^-3]
    inline constexpr double Nv_300K = 1.04e19;   // [cm^-3]

    // Varshni:  E_g(T) = E_g(0) - a T^2 / (T + b)
    inline constexpr double Eg0_Si    = 1.170;   // [eV]
    inline constexpr double varshni_a = 4.73e-4; // [eV/K]
    inline constexpr double varshni_b = 636.0;   // [K]

    // Shallow dopant ionization offsets from the nearest band edge.
    inline constexpr double E_donor_offset    = 0.045;   // P  in Si [eV]
    inline constexpr double E_acceptor_offset = 0.045;   // B  in Si [eV]

    // Degeneracy factors for the ionization statistics.
    inline constexpr double g_donor    = 2.0;
    inline constexpr double g_acceptor = 4.0;

    // ---- Matthiessen mobility parameters (Phase 3 default model) -----------
    //
    //   mu_L(T)   = mu_L(300) * (T / 300)^(-3/2)
    //   mu_I(T,N) = K_I * (T / 300)^( 3/2) * (N_ref / N)
    //   1/mu      = 1/mu_L + 1/mu_I
    //
    // Constants below were tuned so that at T=300 K the model reproduces
    // the textbook Si mobilities
    //     mu_n(pure Si)          ~= 1414 cm^2/Vs
    //     mu_n(N_d = 1e17 cm^-3) ~=  800 cm^2/Vs
    //     mu_p(pure Si)          ~=  470 cm^2/Vs
    //     mu_p(N_a = 1e17 cm^-3) ~=  320 cm^2/Vs
    inline constexpr double mu_L_n_300 = 1414.0;   // lattice-limited electron  mu
    inline constexpr double mu_L_p_300 =  470.5;   // lattice-limited hole      mu
    inline constexpr double mu_I_n_300 = 1786.0;   // impurity-limited electron mu at N=N_ref
    inline constexpr double mu_I_p_300 = 1003.0;   // impurity-limited hole     mu at N=N_ref
    inline constexpr double N_ref_matt = 1.0e17;   // reference impurity conc.  [cm^-3]

    // ---- Arora mobility parameters (Phase 2 legacy, still exposed) ---------
    inline constexpr double mu_n_min_300 = 88.0;
    inline constexpr double mu_n_max_300 = 1414.0;
    inline constexpr double Nref_n_300   = 1.26e17;
    inline constexpr double alpha_n_300  = 0.88;
    inline constexpr double mu_p_min_300 = 54.3;
    inline constexpr double mu_p_max_300 = 470.5;
    inline constexpr double Nref_p_300   = 2.35e17;
    inline constexpr double alpha_p_300  = 0.88;

    // ---- Optical-generation model ------------------------------------------
    //
    //   Above-bandgap illumination creates electron-hole pairs at a volumetric
    //   rate G_opt = Phi * alpha * eta. In steady state the minority-carrier
    //   continuity equation gives
    //
    //       Delta n = Delta p = G_opt * tau
    //
    //   where tau is the recombination lifetime. Instead of tracking Phi,
    //   alpha and tau separately we bundle their product into a single
    //   pedagogical constant K_opt with units of (cm^-3 per eV of excess
    //   photon energy above the gap). The resulting Delta n is then
    //
    //       Delta n = K_opt * (E_photon - E_g)     (hv > E_g)
    //       Delta n = 0                            (hv <= E_g)
    //
    //   K_opt was chosen so that a typical visible-light photon (hv ~ 2.5 eV)
    //   on Si produces Delta n ~ 1e15 cm^-3, consistent with 1-sun injection
    //   levels in lightly doped material.
    inline constexpr double K_opt_excess = 1.0e15;    // [cm^-3 / eV]

    // Clamp wavelength input to a physically sensible window.
    inline constexpr double lambda_min_nm =  300.0;   // near-UV
    inline constexpr double lambda_max_nm = 1500.0;   // short-wave IR

    // Magnetic-field range supported by the UI.
    inline constexpr double B_min_T = 0.0;
    inline constexpr double B_max_T = 10.0;

} // namespace phys


// =============================================================================
// Doping regime selected by the UI.
// =============================================================================
enum class DopingType {
    Intrinsic,   // Pure Si
    NType,       // Phosphorus donors    (N_d)
    PType        // Boron acceptors      (N_a)
};


// =============================================================================
// Mobility model selector (Phase 3).
// =============================================================================
enum class MobilityModel {
    Matthiessen,   // 1/mu = 1/mu_L(T) + 1/mu_I(T,N)
    Arora          // empirical Arora-Hauser-Roulston
};


// =============================================================================
// PhysicsEngine
// =============================================================================
class PhysicsEngine {
public:
    PhysicsEngine();

    // ------------------------- Inputs ---------------------------------------
    void setTemperature          (double T_kelvin);
    void setDopingType           (DopingType type);
    void setDopingConcentration  (double N_per_cm3);
    void setIncompleteIonization (bool enabled);
    void setMobilityModel        (MobilityModel m);

    // Optical
    void setWavelengthNm         (double lambda_nm);
    void setOpticalEnabled       (bool enabled);

    // Magnetic
    void setMagneticField        (double B_tesla);

    // ------------------------- Query: user inputs ---------------------------
    [[nodiscard]] double        getTemperature()           const noexcept { return m_T;                    }
    [[nodiscard]] DopingType    getDopingType()            const noexcept { return m_dopingType;           }
    [[nodiscard]] double        getDopingConcentration()   const noexcept { return m_N;                    }
    [[nodiscard]] bool          getIncompleteIonization()  const noexcept { return m_incompleteIonization; }
    [[nodiscard]] MobilityModel getMobilityModel()         const noexcept { return m_mobilityModel;        }
    [[nodiscard]] double        getWavelengthNm()          const noexcept { return m_lambda_nm;            }
    [[nodiscard]] bool          getOpticalEnabled()        const noexcept { return m_opticalEnabled;       }
    [[nodiscard]] double        getMagneticField()         const noexcept { return m_B;                    }

    // ------------------------- Query: thermodynamic state -------------------
    [[nodiscard]] double getBandgap()               const noexcept { return m_Eg; }
    [[nodiscard]] double getEffectiveNc()           const noexcept { return m_Nc; }
    [[nodiscard]] double getEffectiveNv()           const noexcept { return m_Nv; }
    [[nodiscard]] double getIntrinsicCarrier()      const noexcept { return m_ni; }
    [[nodiscard]] double getElectronConcentration() const noexcept { return m_n;  }
    [[nodiscard]] double getHoleConcentration()     const noexcept { return m_p;  }
    [[nodiscard]] double getFermiLevel()            const noexcept { return m_Ef; }

    // ------------------------- Query: ionization ----------------------------
    [[nodiscard]] double getIonizedDonors()      const noexcept { return m_NdPlus;  }
    [[nodiscard]] double getIonizedAcceptors()   const noexcept { return m_NaMinus; }
    [[nodiscard]] double getIonizationFraction() const noexcept;

    // ------------------------- Query: transport -----------------------------
    [[nodiscard]] double getElectronMobility() const noexcept { return m_mu_n;  }  // [cm^2/Vs]
    [[nodiscard]] double getHoleMobility()     const noexcept { return m_mu_p;  }  // [cm^2/Vs]
    [[nodiscard]] double getConductivity()     const noexcept { return m_sigma; }  // [S/cm]
    [[nodiscard]] double getResistivity()      const noexcept;                     // [Ohm cm]

    // Quasi-total carrier densities including optical excess.
    [[nodiscard]] double getTotalElectronConc() const noexcept { return m_n + m_deltaN; }
    [[nodiscard]] double getTotalHoleConc()     const noexcept { return m_p + m_deltaN; }

    // ------------------------- Query: optical -------------------------------
    [[nodiscard]] double getPhotonEnergy()         const noexcept { return m_Ephoton; } // [eV]
    [[nodiscard]] double getOpticalGeneration()    const noexcept { return m_Gopt;    } // [cm^-3/s]
    [[nodiscard]] double getExcessCarrierDensity() const noexcept { return m_deltaN;  } // [cm^-3]
    [[nodiscard]] bool   isOpticallyPumped()       const noexcept;

    // ------------------------- Query: Hall / magnetic -----------------------
    [[nodiscard]] double getHallCoefficient() const noexcept { return m_R_H; }      // [cm^3/C]
    [[nodiscard]] double getHallVoltage(double current_A,
                                        double thickness_cm,
                                        double field_T) const noexcept;            // [V]

    // ------------------------- Band / dopant levels -------------------------
    [[nodiscard]] double getValenceBandEdge()    const noexcept { return 0.0;  }
    [[nodiscard]] double getConductionBandEdge() const noexcept { return m_Eg; }
    [[nodiscard]] double getDonorLevel()         const noexcept { return m_Eg - phys::E_donor_offset;   }
    [[nodiscard]] double getAcceptorLevel()      const noexcept { return       phys::E_acceptor_offset; }

    // Occupation probability at energy E [eV, relative to E_v].
    [[nodiscard]] double fermiDirac(double E) const noexcept;

    // ------------------------------------------------------------------------
    // CSV export: appends a single row describing the current state.
    // Writes the header row the first time the file is created.
    // ------------------------------------------------------------------------
    [[nodiscard]] bool exportCSV(const std::string& path) const;

    // ------------------------------------------------------------------------
    // Static mobility helpers (exposed for unit tests).
    // ------------------------------------------------------------------------
    [[nodiscard]] static double matthiessenMobilityElectron(double T, double N) noexcept;
    [[nodiscard]] static double matthiessenMobilityHole    (double T, double N) noexcept;
    [[nodiscard]] static double aroraMobilityElectron      (double T, double N) noexcept;
    [[nodiscard]] static double aroraMobilityHole          (double T, double N) noexcept;

    // Helper: photon energy in eV for a given wavelength in nm.
    [[nodiscard]] static double photonEnergyEv(double lambda_nm) noexcept;

private:
    // Driver & solvers
    void recompute();
    void solveFullIonization();
    void solveIncompleteIonization();
    void computeOptical();
    void computeTransport();
    void computeHall();

    // ---- Inputs ------------------------------------------------------------
    double        m_T;           // [K]
    DopingType    m_dopingType;
    double        m_N;           // [cm^-3]
    bool          m_incompleteIonization;
    MobilityModel m_mobilityModel;

    double        m_lambda_nm;   // [nm]
    bool          m_opticalEnabled;

    double        m_B;           // [T]

    // ---- Outputs -----------------------------------------------------------
    double m_Eg, m_Nc, m_Nv, m_ni;
    double m_n,  m_p,  m_Ef;
    double m_NdPlus, m_NaMinus;
    double m_mu_n, m_mu_p, m_sigma;

    double m_Ephoton;   // [eV]
    double m_Gopt;      // [cm^-3/s]     -- stored for reference only
    double m_deltaN;    // [cm^-3]       -- steady-state optical excess

    double m_R_H;       // [cm^3/C]
};
