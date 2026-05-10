// =============================================================================
// Material.cpp
//
//   Author : dex / cemfly-april2026
//   License: MIT
// =============================================================================

#include "Material.hpp"


namespace material {

// =============================================================================
// Silicon (indirect bandgap)
// -----------------------------------------------------------------------------
//   * Workhorse semiconductor for digital electronics.
//   * Indirect gap (Gamma -> X) requires phonon assistance for absorption,
//     so the optical absorption onset is gradual.
// =============================================================================
const Profile Silicon{
    .name              = "Si",
    .donorSpecies      = "P",
    .acceptorSpecies   = "B",
    .isDirectBandgap   = false,

    .Eg0               = 1.170,      // eV at T = 0
    .varshni_a         = 4.73e-4,    // eV/K
    .varshni_b         = 636.0,      // K

    .Nc_300K           = 2.80e19,
    .Nv_300K           = 1.04e19,

    .E_donor_offset    = 0.045,      // P in Si
    .E_acceptor_offset = 0.045,      // B in Si

    .mu_L_n_300        = 1414.0,
    .mu_L_p_300        =  470.5,
    .mu_I_n_300        = 1786.0,
    .mu_I_p_300        = 1003.0,
    .N_ref_matt        = 1.0e17,

    .K_opt_excess      = 1.0e15,

    // High-field saturation (Caughey-Thomas, tabulated in Sze 1.5.4).
    .v_sat_n           = 1.07e7,     // cm/s
    .v_sat_p           = 8.34e6,
    .beta_n            = 2.0,
    .beta_p            = 1.0,

    // SRH lifetimes -- typical bulk float-zone Si (Pierret Ch. 5).
    .tau_n             = 1.0e-3,     // s
    .tau_p             = 1.0e-4,     // s (lower for holes due to common deep-level traps)

    // Early voltage typical for a textbook NPN Si BJT.
    .V_Early           = 75.0,       // V

    .kappa             = 1.50,       // W/(cm K)
    .rho_cp            = 1.65,       // J/(cm^3 K)

    .atomR = 120, .atomG = 170, .atomB = 220,
};

// =============================================================================
// Gallium Arsenide (direct bandgap)
// -----------------------------------------------------------------------------
//   * Direct bandgap at the Gamma point -> efficient radiative recombination
//     (used for laser diodes, LEDs, solar cells) and *much* steeper optical
//     absorption than Si just above E_g.
//   * Electron mobility is ~6x Si's; great for high-frequency transistors.
//   * Asymmetric DOS (Nc << Nv): light electrons, heavy holes.
// =============================================================================
const Profile GalliumArsenide{
    .name              = "GaAs",
    .donorSpecies      = "Si",       // Si on Ga site is a shallow donor
    .acceptorSpecies   = "Be",       // Be is a common shallow acceptor
    .isDirectBandgap   = true,

    .Eg0               = 1.519,      // eV at T = 0
    .varshni_a         = 5.405e-4,
    .varshni_b         = 204.0,

    .Nc_300K           = 4.70e17,
    .Nv_300K           = 9.00e18,

    .E_donor_offset    = 0.006,      // Si donor very shallow in GaAs
    .E_acceptor_offset = 0.028,      // Be acceptor

    .mu_L_n_300        = 8500.0,     // very high electron mobility
    .mu_L_p_300        =  400.0,
    .mu_I_n_300        = 8000.0,
    .mu_I_p_300        =  350.0,
    .N_ref_matt        = 1.0e17,

    // Direct-gap absorption is much stronger -> bigger K_opt.
    .K_opt_excess      = 1.0e16,

    // GaAs displays negative differential mobility above ~3 kV/cm (Gunn);
    // for the basic Caughey-Thomas form we still use a saturation velocity.
    .v_sat_n           = 7.7e6,
    .v_sat_p           = 9.0e6,
    .beta_n            = 2.0,
    .beta_p            = 1.0,

    // Direct-gap GaAs is dominated by radiative recombination, which is
    // captured here with a much shorter SRH lifetime than Si.
    .tau_n             = 5.0e-9,     // s  (~ns regime)
    .tau_p             = 5.0e-9,

    // GaAs HBTs show very high V_A; pedagogical mid-range value.
    .V_Early           = 80.0,

    .kappa             = 0.55,       // GaAs: ~3x lower than Si
    .rho_cp            = 1.74,

    .atomR = 200, .atomG = 130, .atomB = 240,
};

// =============================================================================
// Germanium (indirect bandgap, narrow gap)
// -----------------------------------------------------------------------------
//   * Lower bandgap than Si -> sensitive to longer wavelengths (used for
//     near-IR detectors).
//   * Higher hole mobility than Si.
// =============================================================================
const Profile Germanium{
    .name              = "Ge",
    .donorSpecies      = "P",
    .acceptorSpecies   = "B",
    .isDirectBandgap   = false,

    .Eg0               = 0.7437,     // eV at T = 0
    .varshni_a         = 4.774e-4,
    .varshni_b         = 235.0,

    .Nc_300K           = 1.04e19,
    .Nv_300K           = 6.00e18,

    .E_donor_offset    = 0.012,
    .E_acceptor_offset = 0.011,

    .mu_L_n_300        = 3900.0,
    .mu_L_p_300        = 1900.0,
    .mu_I_n_300        = 4500.0,
    .mu_I_p_300        = 2100.0,
    .N_ref_matt        = 1.0e17,

    .K_opt_excess      = 5.0e15,

    .v_sat_n           = 6.0e6,
    .v_sat_p           = 6.0e6,
    .beta_n            = 2.0,
    .beta_p            = 1.0,

    // Ge has a notably long minority-carrier lifetime in pure crystals.
    .tau_n             = 5.0e-4,
    .tau_p             = 5.0e-4,

    .V_Early           = 50.0,

    .kappa             = 0.60,
    .rho_cp            = 1.65,

    .atomR = 230, .atomG = 160, .atomB =  90,
};

// =============================================================================
// Lookup
// =============================================================================
const Profile& byKind(Kind k) noexcept {
    switch (k) {
        case Kind::Silicon:    return Silicon;
        case Kind::GaAs:       return GalliumArsenide;
        case Kind::Germanium:  return Germanium;
    }
    return Silicon; // unreachable; appeases compiler
}

} // namespace material
