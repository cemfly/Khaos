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

} // namespace material
