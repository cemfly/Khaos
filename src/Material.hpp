#pragma once

// =============================================================================
// Material.hpp -- Si / GaAs / Ge profiles
//
//   Author : dex / cemfly-april2026
//   License: MIT
// -----------------------------------------------------------------------------
// Bundles every material-dependent constant the PhysicsEngine consumes into
// a single immutable struct. Three profiles ship out-of-the-box:
//
//    * Silicon (Si)            -- indirect bandgap, the workhorse
//    * Gallium Arsenide (GaAs) -- direct bandgap (laser/LED material), high
//                                 electron mobility, strong above-gap optical
//                                 absorption
//    * Germanium (Ge)          -- indirect bandgap, narrow gap, used as a
//                                 reference for IR detectors
//
// Sources: Sze "Physics of Semiconductor Devices" (3rd ed.), Kasap (4th ed.),
//          Ioffe Physical-Technical Institute Si/GaAs/Ge handbook.
// =============================================================================

#include <concepts>
#include <string_view>


namespace material {

// -----------------------------------------------------------------------------
// All numeric parameters needed by PhysicsEngine, packed into one POD.
// Units are kept consistent with the rest of the simulator:
//   energies in eV, concentrations in cm^-3, mobilities in cm^2/(V s).
// -----------------------------------------------------------------------------
struct Profile {
    std::string_view name;          // human-readable label, e.g. "Si"
    std::string_view donorSpecies;  // e.g. "P"
    std::string_view acceptorSpecies; // e.g. "B"
    bool   isDirectBandgap;         // true = direct (GaAs), false = indirect

    // ---- Bandgap E_g(T)  =  Eg0 - a T^2 / (T + b)   (Varshni equation) ----
    double Eg0;          // [eV] at T = 0
    double varshni_a;    // [eV/K]
    double varshni_b;    // [K]

    // ---- Effective DOS at 300 K (scales as T^(3/2)) -----------------------
    double Nc_300K;      // [cm^-3]
    double Nv_300K;      // [cm^-3]

    // ---- Shallow dopant ionization energies (offset from band edge) --------
    double E_donor_offset;     // [eV]   E_d = E_c - this
    double E_acceptor_offset;  // [eV]   E_a = E_v + this

    // ---- Matthiessen mobility coefficients --------------------------------
    //
    //   mu_L(T)   = mu_L_300 * (T/300)^(-3/2)               (lattice)
    //   mu_I(T,N) = mu_I_300 * (T/300)^( 3/2) * (N_ref/N)   (impurity)
    //   1/mu      = 1/mu_L + 1/mu_I
    //
    double mu_L_n_300;   // electron lattice-limited mobility at 300 K
    double mu_L_p_300;   // hole     lattice-limited mobility at 300 K
    double mu_I_n_300;   // electron impurity mobility at 300 K, N=N_ref
    double mu_I_p_300;   // hole     impurity mobility at 300 K, N=N_ref
    double N_ref_matt;   // reference impurity concentration

    // ---- Optical generation -----------------------------------------------
    //
    //   Above-gap absorption coefficient depends on whether the bandgap is
    //   direct (alpha ~ (hv - Eg)^(1/2)) or indirect (alpha ~ (hv - Eg)^2,
    //   approximated linearly here for stability over the full slider range).
    //
    //   We bundle the photon flux x absorption x lifetime product into a
    //   single constant K_opt so that
    //
    //       Delta n  =  K_opt * (hv - Eg)^p     (p depends on direct/indirect)
    //
    //   K_opt is tuned per material so GaAs's signature "carrier explosion"
    //   under above-gap illumination is visually obvious next to Si.
    double K_opt_excess; // pedagogical pre-factor

    // ---- High-field transport (Caughey-Thomas saturation) -----------------
    //
    //   v(E) = mu_low * E / [1 + (mu_low * E / v_sat)^beta]^(1/beta)
    //
    //   Phenomenological model from Caughey & Thomas, Proc. IEEE 55 (1967);
    //   widely tabulated in Sze, "Physics of Semiconductor Devices" Sec. 1.5.4.
    //   Si: v_sat,n ~= 1.07e7 cm/s (beta_n=2), v_sat,p ~= 8.34e6 (beta_p=1).
    //   GaAs has negative differential mobility (Gunn effect) above ~3 kV/cm;
    //   we keep the same form for pedagogical purposes.
    double v_sat_n;       // [cm/s] electron saturation velocity
    double v_sat_p;       // [cm/s] hole     saturation velocity
    double beta_n;        // Caughey-Thomas exponent (electrons)
    double beta_p;        // Caughey-Thomas exponent (holes)

    // ---- SRH recombination lifetimes --------------------------------------
    //
    //   R_SRH = (np - n_i^2) / [tau_p (n + n_t) + tau_n (p + p_t)]
    //
    //   with n_t = n_i exp((E_t - E_i)/kT), p_t = n_i exp((E_i - E_t)/kT).
    //   For midgap traps (E_t ~= E_i): n_t = p_t = n_i, simplifying to
    //   R = (np - n_i^2) / [tau (n + p + 2 n_i)].
    //
    //   Pierret Ch. 5 / Sze Sec. 1.5.5. Lifetimes vary by orders of
    //   magnitude across materials and crystal quality.
    double tau_n;         // [s] electron SRH lifetime (typical bulk)
    double tau_p;         // [s] hole     SRH lifetime

    // ---- BJT-relevant (textbook NPN) --------------------------------------
    //
    //   Early voltage governs the slope of I_C vs V_CE in the active region:
    //   I_C(V_CE) = I_C0 * (1 + V_CE / V_A).  A large V_A means a flat
    //   output characteristic (good current source); a small V_A signals
    //   strong base-width modulation. Sze Sec. 5.2.
    double V_Early;       // [V] typical Early voltage for an NPN BJT

    // ---- Thermal (2D heat-equation solver) --------------------------------
    //
    //   Fourier's law of heat conduction:
    //
    //       rho Cp  d T / d t  =  div(kappa grad T) + H_gen(x, y)
    //
    //   On the visualisation grid (dimensionless cells), this collapses to
    //   an FTCS update governed by the diffusivity  alpha = kappa / (rho Cp).
    //
    //   Values are intentionally given in CGS-cm units for consistency with
    //   the rest of the engine (cm, eV, S/cm, ...). The numerical solver
    //   rescales them onto a per-frame coefficient internally.
    //
    //   Si:    kappa ~ 1.5  W/(cm K),  rho_cp ~ 1.65 J/(cm^3 K)
    //   GaAs:  kappa ~ 0.55,           rho_cp ~ 1.74
    //   Ge:    kappa ~ 0.6,            rho_cp ~ 1.65
    double kappa;        // [W / (cm K)]   thermal conductivity at 300 K
    double rho_cp;       // [J / (cm^3 K)] volumetric heat capacity

    // ---- Electrostatic (Poisson solver) ----------------------------------
    //
    //   Static dielectric constant. Required by the Poisson equation
    //
    //       div(eps_s grad psi) = -q (p - n + Nd+ - Na-)
    //
    //   Si:   11.7    GaAs: 12.9    Ge: 16.0    (Sze Tab. 1.1)
    double epsilon_r;    // [-] static relative permittivity

    // ---- Electron affinity (Anderson rule, heterojunctions)  [Phase 5] --
    //
    //   chi  = energy required to lift a conduction-band electron to the
    //          vacuum level [eV].  Together with E_g it pins the band
    //          edges in an absolute reference frame:
    //
    //              E_c(x) = E_vac(x) - chi(x)  =  -q psi(x) - chi(x)
    //              E_v(x) = E_c(x) - E_g(x)
    //
    //   At an A/B heterojunction the band edges step by:
    //
    //              Delta E_c =  chi_A - chi_B
    //              Delta E_v = (chi_B - chi_A) + (E_g,B - E_g,A)
    //                        = -Delta E_c + Delta E_g
    //
    //   Tabulated (300 K, Sze Sec. 1.5 / Pierret App. F):
    //     Si   :  chi = 4.05 eV    Ge   :  chi = 4.00 eV
    //     GaAs :  chi = 4.07 eV
    //
    //   With these three numbers the engine reproduces the textbook
    //   Si/Ge type-II offset (Delta E_v ~ 0.42 eV, Delta E_c ~ 0.05 eV)
    //   and the Si/GaAs broken-gap alignment.
    double chi;          // [eV] electron affinity

    // ---- Auger recombination ---------------------------------------------
    //
    //   R_Aug = (C_n n + C_p p) (n p - n_i^2)
    //
    //   Three-particle process: an excess pair recombines and the released
    //   energy is absorbed by a third carrier (electron-electron-hole or
    //   hole-hole-electron). Dominates SRH/radiative for N > ~1e18 cm^-3
    //   and is the principal loss mechanism in solar cell emitters and
    //   heavily doped BJT bases.  Sze Sec. 1.5.6 / Pierret Sec. 5.2.4.
    //
    //   Tabulated values at 300 K (Si: Dziewior & Schmid 1977):
    //     Si   : C_n ~ 2.8e-31 cm^6/s,  C_p ~ 9.9e-32
    //     GaAs : C_n ~ 1.0e-30,         C_p ~ 1.0e-30   (direct gap;
    //                                                     much smaller in
    //                                                     practice -- here
    //                                                     pedagogically
    //                                                     bumped for
    //                                                     visibility)
    //     Ge   : C_n ~ 8.0e-32,         C_p ~ 2.8e-31
    double C_n_aug;      // [cm^6/s] electron Auger coefficient
    double C_p_aug;      // [cm^6/s] hole     Auger coefficient

    // ---- Impact ionization (Chynoweth) -----------------------------------
    //
    //   alpha(E) = alpha_inf * exp(-(E_crit / E)^m)        (Chynoweth, 1958)
    //
    //   alpha [cm^-1] is the ionization coefficient: number of secondary
    //   electron-hole pairs created per cm of carrier travel. Avalanche
    //   breakdown occurs when the ionization integral
    //
    //     integral_0^W alpha_n exp{int_0^x (alpha_p - alpha_n) dx'} dx --> 1
    //
    //   reaches unity. Sze Sec. 2.4.2 / Pierret Sec. 6.2.3.
    //
    //   Si (Maes/Van Overstraeten):
    //       alpha_inf,n ~ 7.03e5 cm^-1, E_crit,n ~ 1.231e6 V/cm, m_n = 1
    //       alpha_inf,p ~ 1.582e6,      E_crit,p ~ 2.036e6,      m_p = 1
    //   GaAs / Ge values follow the Sze tables.
    double alpha_inf_n;  // [cm^-1] electron Chynoweth pre-factor
    double alpha_inf_p;  // [cm^-1] hole     Chynoweth pre-factor
    double E_crit_n;     // [V/cm]  electron critical field
    double E_crit_p;     // [V/cm]  hole     critical field
    double chyn_m;       // [-]     exponent (1 for Si, ~1 for GaAs/Ge)

    // ---- Band-to-band (Zener) tunneling: Kane model -----------------------
    //
    //   G_BTBT(E) = A_kane * E^P * exp(-B_kane / E)            (Kane 1961)
    //
    //   with P = 2 (direct gap) or 5/2 (indirect, phonon-assisted). The
    //   B parameter encodes the tunneling barrier:
    //
    //     B_kane = (pi sqrt(m*) E_g^(3/2)) / (2 sqrt(2) q hbar)
    //
    //   and is therefore strongly material-dependent (small for narrow
    //   gap / low effective mass).
    //
    //   References: Kane J. Phys. Chem. Solids 12 (1959) 181;
    //               Sze Sec. 4.3 (Zener diodes); Hurkx et al. IEEE-ED 39 (1992).
    //
    //   Tabulated:
    //     Si   : A ~ 3.5e21 (cm^-3 s^-1)/(V/cm)^P,  B ~ 2.25e7 V/cm
    //     GaAs : A ~ 1.0e20 (direct, P=2),          B ~ 4.30e7
    //     Ge   : A ~ 1.6e22,                         B ~ 1.20e7
    bool   btbt_isDirect; // power P = 2 if direct, 5/2 if indirect
    double A_kane;        // [cm^-3 s^-1 (V/cm)^-P]
    double B_kane;        // [V/cm]

    // ---- Display ----------------------------------------------------------
    // Atom rendering colour (Si-blue, GaAs-purple, Ge-orange). Stored as raw
    // RGB bytes so this header has zero SFML dependency.
    unsigned char atomR;
    unsigned char atomG;
    unsigned char atomB;
};


// =============================================================================
// Material catalogue
// =============================================================================
extern const Profile Silicon;
extern const Profile GalliumArsenide;
extern const Profile Germanium;


enum class Kind { Silicon = 0, GaAs = 1, Germanium = 2 };

// Returns one of the three profiles above.
const Profile& byKind(Kind k) noexcept;

// Convenience: stable list for UI dropdowns.
inline constexpr int  kCount        = 3;
inline constexpr const char* kLabels[kCount] = { "Silicon (Si)",
                                                  "Gallium Arsenide (GaAs)",
                                                  "Germanium (Ge)" };


// =============================================================================
// C++20 concept -- compile-time validation of "anything that looks like a
// semiconductor material profile". Static helpers in PhysicsEngine accept
// types that satisfy this concept, so a wrong-type argument fails to compile
// with a clear "constraint not satisfied" diagnostic instead of opaque
// member-access errors.
// =============================================================================
template <typename T>
concept SemiconductorProfile = requires(const T& t) {
    { t.Eg0              } -> std::convertible_to<double>;
    { t.varshni_a        } -> std::convertible_to<double>;
    { t.varshni_b        } -> std::convertible_to<double>;
    { t.Nc_300K          } -> std::convertible_to<double>;
    { t.Nv_300K          } -> std::convertible_to<double>;
    { t.isDirectBandgap  } -> std::convertible_to<bool>;
    { t.v_sat_n          } -> std::convertible_to<double>;
    { t.tau_n            } -> std::convertible_to<double>;
    { t.epsilon_r        } -> std::convertible_to<double>;
    { t.chi              } -> std::convertible_to<double>;
    { t.C_n_aug          } -> std::convertible_to<double>;
    { t.alpha_inf_n      } -> std::convertible_to<double>;
    { t.E_crit_n         } -> std::convertible_to<double>;
    { t.A_kane           } -> std::convertible_to<double>;
    { t.B_kane           } -> std::convertible_to<double>;
};

// Sanity: Profile itself must satisfy the concept it advertises.
static_assert(SemiconductorProfile<Profile>);

} // namespace material
