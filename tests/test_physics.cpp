// =============================================================================
// tests/test_physics.cpp
//
//   Author : dex / cemfly-april2026
//   License: MIT
// -----------------------------------------------------------------------------
// GoogleTest-based unit tests for PhysicsEngine.
//
// Covers every physics regime shipped with the simulator:
//
//   Thermal
//     - Intrinsic carrier concentration at 300 K
//     - Mass action law  n p = n_i^2
//     - Varshni bandgap monotonically decreases with T
//
//   Ionization
//     - Shallow donors fully ionised at 300 K
//     - Freeze-out at cryogenic T (50 K)
//
//   Transport
//     - Matthiessen mobility reproduces textbook room-T values
//     - Matthiessen mobility drops with doping (impurity scattering dominates)
//     - Matthiessen lattice limit:  mu propto T^(-3/2)
//     - Arora model still matches textbook values (regression test)
//     - Conductivity increases with doping and with illumination
//
//   Optical
//     - Photon energy obeys  E = h c / lambda
//     - Delta n = 0 when hv <= E_g (no absorption)
//     - Delta n > 0 when hv >  E_g (pair generation)
//
//   Magnetic / Hall
//     - R_H < 0 for n-type; R_H > 0 for p-type
//     - |R_H| scales as 1/(n q) in the single-carrier limit
//
// Build:
//     cmake -S . -B build -DBUILD_TESTS=ON
//     cmake --build build --target physics_tests
//     ctest --test-dir build --output-on-failure
// =============================================================================

#include <gtest/gtest.h>

#include <cmath>

#include "PhysicsEngine.hpp"


// =============================================================================
// Thermal
// =============================================================================

TEST(Thermal, IntrinsicCarrierAt300K) {
    PhysicsEngine pe;
    pe.setTemperature(300.0);
    pe.setDopingType(DopingType::Intrinsic);

    const double ni = pe.getIntrinsicCarrier();
    EXPECT_GT(ni, 1.0e9)  << "n_i = " << ni;
    EXPECT_LT(ni, 1.0e11) << "n_i = " << ni;
}


TEST(Thermal, MassActionLaw) {
    PhysicsEngine pe;
    pe.setTemperature(300.0);
    pe.setDopingType(DopingType::NType);
    pe.setDopingConcentration(1.0e16);

    const double n  = pe.getElectronConcentration();
    const double p  = pe.getHoleConcentration();
    const double ni = pe.getIntrinsicCarrier();
    EXPECT_NEAR((n * p) / (ni * ni), 1.0, 1.0e-3);
}


TEST(Thermal, BandgapDecreasesWithTemperature) {
    PhysicsEngine pe;
    pe.setTemperature(100.0);
    const double Eg_low = pe.getBandgap();
    pe.setTemperature(500.0);
    const double Eg_high = pe.getBandgap();

    EXPECT_GT(Eg_low, Eg_high);
    EXPECT_GT(Eg_low,  1.05);
    EXPECT_LT(Eg_high, 1.20);
}


// =============================================================================
// Ionization
// =============================================================================

TEST(Ionization, FullIonizationAt300K) {
    PhysicsEngine pe;
    pe.setIncompleteIonization(true);
    pe.setDopingType(DopingType::NType);
    pe.setDopingConcentration(1.0e16);
    pe.setTemperature(300.0);

    EXPECT_GT(pe.getIonizationFraction(), 0.95);
}


TEST(Ionization, FreezeOutAt50K) {
    PhysicsEngine pe;
    pe.setIncompleteIonization(true);
    pe.setDopingType(DopingType::NType);
    pe.setDopingConcentration(1.0e16);
    pe.setTemperature(50.0);

    const double n  = pe.getElectronConcentration();
    const double Nd = pe.getDopingConcentration();

    EXPECT_LT(n, 0.5 * Nd);
    EXPECT_LT(pe.getIonizationFraction(), 0.5);

    // Fermi level pinned near the donor level at freeze-out.
    EXPECT_NEAR(pe.getFermiLevel(), pe.getDonorLevel(), 0.05);

    // Sanity: turning freeze-out OFF recovers n ~ N_d.
    pe.setIncompleteIonization(false);
    EXPECT_NEAR(pe.getElectronConcentration() / Nd, 1.0, 0.01);
}


// =============================================================================
// Transport -- Matthiessen model
// =============================================================================

TEST(Mobility, MatthiessenAt300K) {
    const auto& Si = material::Silicon;
    EXPECT_NEAR(PhysicsEngine::matthiessenMobilityElectron(Si, 300.0, 0.0),
                1414.0, 1.0);
    EXPECT_NEAR(PhysicsEngine::matthiessenMobilityHole    (Si, 300.0, 0.0),
                 470.5, 1.0);

    const double mu_n_heavy =
        PhysicsEngine::matthiessenMobilityElectron(Si, 300.0, 1.0e19);
    EXPECT_LT(mu_n_heavy, 300.0);
}


TEST(Mobility, MatthiessenTemperatureScalingIntrinsic) {
    const auto& Si = material::Silicon;
    const double mu_300 = PhysicsEngine::matthiessenMobilityElectron(Si, 300.0, 0.0);
    const double mu_600 = PhysicsEngine::matthiessenMobilityElectron(Si, 600.0, 0.0);
    const double expectedRatio = std::pow(600.0 / 300.0, -1.5);
    EXPECT_NEAR(mu_600 / mu_300, expectedRatio, 1.0e-3);
}


TEST(Material, GaAsHasDirectBandgapAndHigherMobility) {
    EXPECT_TRUE (material::GalliumArsenide.isDirectBandgap);
    EXPECT_FALSE(material::Silicon       .isDirectBandgap);
    EXPECT_FALSE(material::Germanium     .isDirectBandgap);

    // GaAs electron mobility is famously much higher than Si.
    EXPECT_GT(material::GalliumArsenide.mu_L_n_300,
              material::Silicon.mu_L_n_300 * 3.0);
}


TEST(Material, OpticalAbsorptionStrongerInGaAs) {
    PhysicsEngine si  (material::Kind::Silicon);
    PhysicsEngine gaas(material::Kind::GaAs);
    for (auto* pe : { &si, &gaas }) {
        pe->setTemperature(300.0);
        pe->setDopingType(DopingType::Intrinsic);
        pe->setOpticalEnabled(true);
        pe->setWavelengthNm(500.0);   // 2.48 eV, well above both gaps
    }
    // Direct GaAs should produce a much larger excess than indirect Si
    // for the same photon energy.
    EXPECT_GT(gaas.getExcessCarrierDensity(),
              5.0 * si.getExcessCarrierDensity());
}


// =============================================================================
// Static T-evaluators (electrothermal feedback path)
// =============================================================================

TEST(Static, BandgapMonotonicityAcrossMaterials) {
    // E_g(0) > E_g(300 K) > E_g(600 K) for every shipped material.
    for (const auto* m : { &material::Silicon,
                            &material::GalliumArsenide,
                            &material::Germanium })
    {
        const double a = PhysicsEngine::bandgapAt(*m,   1.0);
        const double b = PhysicsEngine::bandgapAt(*m, 300.0);
        const double c = PhysicsEngine::bandgapAt(*m, 600.0);
        EXPECT_GT(a, b) << m->name;
        EXPECT_GT(b, c) << m->name;
    }
}


TEST(Static, IntrinsicCarrierExplodesWithTemperature) {
    // Doubling T from 300 K to 600 K should bump n_i by several orders of
    // magnitude in Si due to the exp(-Eg/2kT) factor.
    const auto& Si = material::Silicon;
    const double ni_300 = PhysicsEngine::intrinsicCarrierAt(Si, 300.0);
    const double ni_600 = PhysicsEngine::intrinsicCarrierAt(Si, 600.0);
    EXPECT_GT(ni_600 / ni_300, 100.0);
}


// =============================================================================
// Mobility uses ionized impurity count (freeze-out -> higher mobility)
// =============================================================================

TEST(Mobility, FreezeOutIncreasesMobility) {
    // At cryogenic T, incomplete ionization should free electrons from
    // the impurity-scattering picture, so mu_n should be HIGHER with
    // freeze-out enabled than with full ionization (where every donor
    // still scatters).
    PhysicsEngine pe;
    pe.setDopingType(DopingType::NType);
    pe.setDopingConcentration(1.0e17);
    pe.setTemperature(80.0);

    pe.setIncompleteIonization(false);
    const double mu_full = pe.getElectronMobility();

    pe.setIncompleteIonization(true);
    const double mu_freeze = pe.getElectronMobility();

    EXPECT_GT(mu_freeze, mu_full)
        << "mu_full=" << mu_full << "  mu_freeze=" << mu_freeze;
}


// =============================================================================
// Bandgap value sanity (reproduces textbook Si E_g ~ 1.12 eV at 300 K)
// =============================================================================

TEST(Static, SiBandgapMatchesTextbookAt300K) {
    const double Eg = PhysicsEngine::bandgapAt(material::Silicon, 300.0);
    EXPECT_NEAR(Eg, 1.12, 0.01);     // literature: 1.12 eV
}


TEST(Static, GaAsBandgapMatchesTextbookAt300K) {
    const double Eg = PhysicsEngine::bandgapAt(material::GalliumArsenide, 300.0);
    EXPECT_NEAR(Eg, 1.42, 0.02);     // literature: 1.424 eV
}


// =============================================================================
// Numerical stability: log-domain n_i survives cryogenic temperatures
// =============================================================================

TEST(Numerical, LogIntrinsicCarrierStableAt10K) {
    // Linear form underflows  (n_i ~ 10^-280 cm^-3 << denormal range).
    // The log10 form must still return a finite, very negative number.
    const double l10 = PhysicsEngine::log10IntrinsicCarrierAt(
        material::Silicon, 10.0);
    EXPECT_LT(l10, -200.0);
    EXPECT_GT(l10, -400.0);
    EXPECT_FALSE(std::isnan(l10));
    EXPECT_FALSE(std::isinf(l10));
}


TEST(Numerical, IntrinsicCarrierClampsToZeroNotNaN) {
    // Linear evaluator must not return NaN at extreme cryogenic T.
    const double ni = PhysicsEngine::intrinsicCarrierAt(material::Silicon, 5.0);
    EXPECT_FALSE(std::isnan(ni));
    EXPECT_GE(ni, 0.0);
}


// =============================================================================
// Caughey-Thomas high-field saturation
// =============================================================================

TEST(HighField, MobilityRecoversLowFieldLimit) {
    // At E -> 0, mu_eff must equal mu_low.
    const double mu = PhysicsEngine::highFieldMobility(
        1414.0, 0.0, 1.07e7, 2.0);
    EXPECT_NEAR(mu, 1414.0, 1e-9);
}


TEST(HighField, VelocitySaturatesAtHighField) {
    // At E -> infinity, drift velocity saturates at v_sat (cm/s).
    constexpr double mu_low = 1414.0;        // cm^2/Vs
    constexpr double v_sat  = 1.07e7;        // cm/s
    const double E      = 1.0e6;             // V/cm  (huge)
    const double mu_eff = PhysicsEngine::highFieldMobility(
        mu_low, E, v_sat, 2.0);
    const double v      = mu_eff * E;
    EXPECT_NEAR(v, v_sat, 0.05 * v_sat);     // within 5% of saturation
    EXPECT_LT (mu_eff, mu_low);              // strictly degraded
}


// =============================================================================
// SRH recombination
// =============================================================================

TEST(SRH, ZeroAtEquilibrium) {
    // np = n_i^2  =>  R = 0  (detailed balance).
    constexpr double n_i = 1.0e10;
    const double R = PhysicsEngine::recombSRH(
        n_i, n_i, n_i, 1.0e-6, 1.0e-6);
    EXPECT_NEAR(R, 0.0, 1.0e-3);
}


TEST(SRH, PositiveUnderInjection) {
    // np > n_i^2 (e.g. illuminated)  =>  R > 0  (net recombination).
    constexpr double n_i = 1.0e10;
    const double R = PhysicsEngine::recombSRH(
        1.0e15, 1.0e15, n_i, 1.0e-6, 1.0e-6);
    EXPECT_GT(R, 0.0);
}


// =============================================================================
// Bandgap narrowing (Slotboom)
// =============================================================================

TEST(BandgapNarrowing, NegligibleAtLowDoping) {
    // Slotboom's expression is continuous, so it produces a sub-meV tail
    // at N << N_ref instead of an exact zero. "Negligible" relative to a
    // ~1 eV bandgap means < 1 meV, which is the tolerance we use here.
    EXPECT_LT(PhysicsEngine::bandgapNarrowing(1.0e15), 1.0e-3);
}


TEST(BandgapNarrowing, GrowsLogarithmicallyForHeavyDoping) {
    const double d18 = PhysicsEngine::bandgapNarrowing(1.0e18);
    const double d19 = PhysicsEngine::bandgapNarrowing(1.0e19);
    const double d20 = PhysicsEngine::bandgapNarrowing(1.0e20);
    EXPECT_GT(d18, 0.0);
    EXPECT_GT(d19, d18);
    EXPECT_GT(d20, d19);
    EXPECT_LT(d20, 0.20);                    // < 200 meV is the realistic bound
}


// =============================================================================
// Optical absorption + penetration depth
// =============================================================================

TEST(Optical, GaAsAbsorptionMuchStrongerThanSi) {
    constexpr double hv = 2.5;               // ~ 500 nm visible
    const double a_si  = PhysicsEngine::absorptionCoefficient(
        material::Silicon, hv);
    const double a_gaas = PhysicsEngine::absorptionCoefficient(
        material::GalliumArsenide, hv);
    EXPECT_GT(a_gaas, 5.0 * a_si);
}


TEST(Optical, NoAbsorptionBelowGap) {
    EXPECT_EQ(PhysicsEngine::absorptionCoefficient(material::Silicon, 0.5),
              0.0);
    EXPECT_TRUE(std::isinf(
        PhysicsEngine::penetrationDepthCm(material::Silicon, 0.5)));
}


// =============================================================================
// BJT injection efficiency / current gain / Early effect
// =============================================================================

TEST(BJT, EmitterEfficiencyApproachesUnityForHeavyEmitter) {
    // N_E >> N_B  =>  gamma -> 1
    const double gamma = PhysicsEngine::bjtEmitterEfficiency(
        /*D_n=*/35.0, /*N_E=*/1.0e20, /*W_E=*/1.0e-4,
        /*D_p=*/12.0, /*N_B=*/1.0e16, /*W_B=*/0.5e-4);
    EXPECT_GT(gamma, 0.99);
}


TEST(BJT, CurrentGainPositive) {
    // gamma * alpha_T near unity gives large beta.
    const double beta = PhysicsEngine::bjtCurrentGain(0.998, 0.99);
    EXPECT_GT(beta, 80.0);           // textbook range
    EXPECT_LT(beta, 5000.0);
}


TEST(BJT, EarlyEffectIncreasesIc) {
    EXPECT_NEAR(PhysicsEngine::earlyEffectFactor(0.0, 75.0), 1.0, 1.0e-9);
    EXPECT_GT  (PhysicsEngine::earlyEffectFactor(5.0, 75.0), 1.0);
}


TEST(Mobility, AroraRegression) {
    // Kept to detect accidental regressions in the legacy Arora model.
    EXPECT_NEAR(PhysicsEngine::aroraMobilityElectron(300.0, 0.0),
                1414.0, 50.0);
    EXPECT_NEAR(PhysicsEngine::aroraMobilityHole    (300.0, 0.0),
                 470.5, 30.0);
}


// =============================================================================
// Transport -- Conductivity
// =============================================================================

TEST(Conductivity, IncreasesWithDoping) {
    PhysicsEngine a, b, c;
    for (auto* pe : {&a, &b, &c}) {
        pe->setTemperature(300.0);
        pe->setDopingType(DopingType::NType);
    }
    a.setDopingConcentration(1.0e14);
    b.setDopingConcentration(1.0e16);
    c.setDopingConcentration(1.0e18);
    EXPECT_LT(a.getConductivity(), b.getConductivity());
    EXPECT_LT(b.getConductivity(), c.getConductivity());
}


TEST(Conductivity, IlluminationBoostsIntrinsicSigma) {
    PhysicsEngine dark, lit;
    for (auto* pe : {&dark, &lit}) {
        pe->setDopingType(DopingType::Intrinsic);
        pe->setTemperature(300.0);
    }
    // Green light ~ 2.5 eV, well above Si bandgap.
    lit.setWavelengthNm(500.0);
    lit.setOpticalEnabled(true);

    EXPECT_GT(lit.getExcessCarrierDensity(), 0.0);
    EXPECT_GT(lit.getConductivity(), dark.getConductivity() * 10.0)
        << "Illumination should dominate over intrinsic carriers in pure Si";
}


// =============================================================================
// Optical absorption
// =============================================================================

TEST(Optical, PhotonEnergyFormula) {
    // E[eV] = 1239.84 / lambda[nm]
    EXPECT_NEAR(PhysicsEngine::photonEnergyEv(500.0),  2.4797, 1.0e-3);
    EXPECT_NEAR(PhysicsEngine::photonEnergyEv(1240.0), 1.0000, 1.0e-3);
}


TEST(Optical, NoAbsorptionBelowBandgap) {
    PhysicsEngine pe;
    pe.setTemperature(300.0);
    pe.setWavelengthNm(1500.0);     // hv ~ 0.83 eV  < E_g
    pe.setOpticalEnabled(true);

    EXPECT_LT(pe.getPhotonEnergy(), pe.getBandgap());
    EXPECT_EQ(pe.getExcessCarrierDensity(), 0.0);
    EXPECT_FALSE(pe.isOpticallyPumped());
}


TEST(Optical, AbsorptionAboveBandgap) {
    PhysicsEngine pe;
    pe.setTemperature(300.0);
    pe.setWavelengthNm(500.0);     // hv ~ 2.48 eV  >> E_g
    pe.setOpticalEnabled(true);

    EXPECT_GT(pe.getPhotonEnergy(), pe.getBandgap());
    EXPECT_GT(pe.getExcessCarrierDensity(), 0.0);
    EXPECT_TRUE(pe.isOpticallyPumped());
}


// =============================================================================
// Magnetic / Hall
// =============================================================================

TEST(Hall, SignIsNegativeForNType) {
    PhysicsEngine pe;
    pe.setTemperature(300.0);
    pe.setDopingType(DopingType::NType);
    pe.setDopingConcentration(1.0e17);

    EXPECT_LT(pe.getHallCoefficient(), 0.0);
}


TEST(Hall, SignIsPositiveForPType) {
    PhysicsEngine pe;
    pe.setTemperature(300.0);
    pe.setDopingType(DopingType::PType);
    pe.setDopingConcentration(1.0e17);

    EXPECT_GT(pe.getHallCoefficient(), 0.0);
}


TEST(Hall, SingleCarrierLimit) {
    // For heavily doped n-type, R_H ~ -1 / (n q).
    PhysicsEngine pe;
    pe.setTemperature(300.0);
    pe.setDopingType(DopingType::NType);
    pe.setDopingConcentration(1.0e18);

    const double n    = pe.getElectronConcentration();
    const double Rh   = pe.getHallCoefficient();
    const double Rh_0 = -1.0 / (n * phys::q_e);    // simple one-band formula

    // Allow a 30% tolerance -- the two-band formula always reduces
    // |R_H| slightly vs. the simple approximation.
    EXPECT_NEAR(Rh / Rh_0, 1.0, 0.30);
}


// =============================================================================
// Fermi level response to doping
// =============================================================================

TEST(FermiLevel, TrackingDoping) {
    PhysicsEngine pe;
    pe.setTemperature(300.0);

    pe.setDopingType(DopingType::NType);
    pe.setDopingConcentration(1.0e17);
    EXPECT_GT(pe.getFermiLevel(), 0.5 * pe.getBandgap());

    pe.setDopingType(DopingType::PType);
    pe.setDopingConcentration(1.0e17);
    EXPECT_LT(pe.getFermiLevel(), 0.5 * pe.getBandgap());
}


// =============================================================================
// Scharfetter-Gummel Bernoulli function  (Phase 4)
// -----------------------------------------------------------------------------
//   B(x) = x / (exp(x) - 1).
//   Identity in equilibrium drift-diffusion: B(0) = 1, B(x) - B(-x) = -x.
//   Asymptotics: B(x)  -> 0   as x -> +inf
//                B(x)  -> -x  as x -> -inf
// =============================================================================
TEST(Bernoulli, IdentityAtZero) {
    EXPECT_NEAR(PhysicsEngine::bernoulli(0.0), 1.0, 1.0e-12);
}

TEST(Bernoulli, IdentityTinyArgs) {
    // Taylor branch (|x| < 1e-6)
    for (double x : {-1e-8, -1e-10, 1e-12, 1e-7}) {
        const double b = PhysicsEngine::bernoulli(x);
        EXPECT_NEAR(b, 1.0 - 0.5 * x, 1.0e-9);
    }
}

TEST(Bernoulli, AntisymmetryIdentity) {
    // For all x != 0:  B(x) - B(-x) = -x   (SG drift-diffusion identity)
    for (double x : {-3.0, -1.0, -0.1, 0.1, 1.0, 3.0, 12.0}) {
        const double lhs = PhysicsEngine::bernoulli( x)
                         - PhysicsEngine::bernoulli(-x);
        EXPECT_NEAR(lhs, -x, 1.0e-9) << "x = " << x;
    }
}

TEST(Bernoulli, OverflowBranchStableAtX50) {
    // x = 50 -> exp(50) ~ 5.18e21; naive (exp(x) - 1) overflows in
    // single-precision but we're in double-precision land. Still: the
    // |x|>40 branch uses asymptotic x*exp(-x) which is finite & small.
    const double b = PhysicsEngine::bernoulli(50.0);
    EXPECT_TRUE(std::isfinite(b));
    EXPECT_GT(b, 0.0);
    EXPECT_LT(b, 1.0e-18);
}

TEST(Bernoulli, UnderflowBranchAtXMinus50) {
    // x = -50: exp(x) ~ 0, denominator ~ -1, so B(x) ~ -x = 50.
    const double b = PhysicsEngine::bernoulli(-50.0);
    EXPECT_NEAR(b, 50.0, 1.0e-9);
}


// =============================================================================
// Recombination triplet  (Phase 2-4)
// -----------------------------------------------------------------------------
// Detailed-balance: every mechanism's rate must vanish when n p = n_i^2.
// Three-particle Auger > SRH at heavy doping; radiative dominates in
// direct-gap GaAs at moderate forward injection.
// Sze 1.5.6 / Pankove Ch. 6.
// =============================================================================
TEST(Recombination, AugerZeroAtEquilibrium) {
    const double n_i = 1.0e10;
    const double R = PhysicsEngine::recombAuger(n_i, n_i, n_i,
                                                2.8e-31, 9.9e-32);
    EXPECT_NEAR(R, 0.0, 1.0e-6);
}

TEST(Recombination, RadiativeZeroAtEquilibrium) {
    const double n_i = 1.0e10;
    const double R = PhysicsEngine::recombRadiative(n_i, n_i, n_i, 1.1e-14);
    EXPECT_NEAR(R, 0.0, 1.0e-6);
}

TEST(Recombination, RadiativeZeroWhenCoefficientIsZero) {
    EXPECT_DOUBLE_EQ(
        PhysicsEngine::recombRadiative(1.0e18, 1.0e18, 1.0e10, 0.0), 0.0);
}

TEST(Recombination, NetIsSumOfThree) {
    // Net = SRH + Auger + Radiative -- verify the aggregator agrees
    // with explicit summation, using the GaAs profile (largest B_rad
    // so radiative contribution is non-trivial).
    const auto& mat = material::byKind(material::Kind::GaAs);
    const double n_i = 2.0e6, n = 1.0e16, p = 1.0e16;
    const double R_srh = PhysicsEngine::recombSRH(n, p, n_i,
                                                  mat.tau_n, mat.tau_p);
    const double R_aug = PhysicsEngine::recombAuger(n, p, n_i,
                                                    mat.C_n_aug, mat.C_p_aug);
    const double R_rad = PhysicsEngine::recombRadiative(n, p, n_i, mat.B_rad);
    const double U     = PhysicsEngine::netRecombination(n, p, n_i, mat);
    EXPECT_NEAR(U, R_srh + R_aug + R_rad,
                std::max(1.0e-6, 1.0e-9 * std::abs(U)));
}

TEST(Recombination, GaAsRadiativeMuchLargerThanSi) {
    // Physical signature of direct vs indirect gap: with the same n,p,
    // GaAs R_rad >> Si R_rad by ~5 orders.
    const double n = 1.0e17, p = 1.0e17, n_i = 1.0e10;
    const double R_Si   = PhysicsEngine::recombRadiative(
        n, p, n_i, material::Silicon.B_rad);
    const double R_GaAs = PhysicsEngine::recombRadiative(
        n, p, 2.0e6, material::GalliumArsenide.B_rad);
    EXPECT_GT(R_GaAs / R_Si, 1.0e3);
}


// =============================================================================
// Impact ionisation (Chynoweth) + band-to-band tunnelling (Kane)  (Phase 2)
// -----------------------------------------------------------------------------
// alpha(E) -> 0 at E = 0, monotonically increasing in |E|.
// G_BTBT  -> 0 at E = 0, also monotonic. Both must be finite and >= 0.
// =============================================================================
TEST(Chynoweth, ZeroAtZeroField) {
    const double a = PhysicsEngine::chynowethRate(
        0.0, 7.03e5, 1.231e6, 1.0);
    EXPECT_DOUBLE_EQ(a, 0.0);
}

TEST(Chynoweth, MonotonicInField) {
    const double a1 = PhysicsEngine::chynowethRate(
        1.0e5, 7.03e5, 1.231e6, 1.0);
    const double a2 = PhysicsEngine::chynowethRate(
        5.0e5, 7.03e5, 1.231e6, 1.0);
    const double a3 = PhysicsEngine::chynowethRate(
        2.0e6, 7.03e5, 1.231e6, 1.0);
    EXPECT_LT(a1, a2);
    EXPECT_LT(a2, a3);
}

TEST(Chynoweth, NonNegativeAcrossRange) {
    for (double E = 0.0; E < 1.0e7; E += 5.0e5) {
        EXPECT_GE(PhysicsEngine::chynowethRate(
            E, 7.03e5, 1.231e6, 1.0), 0.0);
    }
}

TEST(KaneBTBT, ZeroAtZeroField) {
    EXPECT_DOUBLE_EQ(
        PhysicsEngine::kaneBTBT(0.0, 3.5e21, 2.25e7, false), 0.0);
}

TEST(KaneBTBT, IndirectVsDirectAtHighField) {
    // Direct (P=2) and indirect (P=5/2) both monotonic; for the same
    // A,B the indirect form has the extra sqrt(E) factor.
    const double E = 1.5e6;
    const double Gd = PhysicsEngine::kaneBTBT(E, 1.0e21, 2.0e7, true);
    const double Gi = PhysicsEngine::kaneBTBT(E, 1.0e21, 2.0e7, false);
    EXPECT_GT(Gi, Gd);  // sqrt(E) > 1 for E > 1
}

TEST(KaneBTBT, MonotonicInField) {
    using P = PhysicsEngine;
    EXPECT_LT(P::kaneBTBT(1.0e6, 3.5e21, 2.25e7, false),
              P::kaneBTBT(2.0e6, 3.5e21, 2.25e7, false));
    EXPECT_LT(P::kaneBTBT(2.0e6, 3.5e21, 2.25e7, false),
              P::kaneBTBT(5.0e6, 3.5e21, 2.25e7, false));
}


// =============================================================================
// Poisson / electrostatics helpers  (Phase 1, 3)
// -----------------------------------------------------------------------------
//   V_T = k_B T (eV/K * K -> V via numerical equivalence)
//   V_bi = V_T ln(Nd Na / n_i^2)        (Sze Eq. 2.10)
// =============================================================================
TEST(Electrostatics, ThermalVoltageAt300K) {
    const double V_T = PhysicsEngine::thermalVoltage(300.0);
    EXPECT_NEAR(V_T, 0.02585, 1.0e-4);
}

TEST(Electrostatics, BuiltInPotentialFormula) {
    const double V_T  = PhysicsEngine::thermalVoltage(300.0);
    const double n_i  = 1.0e10;
    const double Nd   = 1.0e17;
    const double Na   = 1.0e17;
    const double Vbi  = PhysicsEngine::builtInPotential(Nd, Na, n_i, V_T);
    // V_bi = V_T ln(Nd Na / n_i^2)
    const double expected = V_T * std::log((Nd * Na) / (n_i * n_i));
    EXPECT_NEAR(Vbi, expected, 1.0e-9);
}

TEST(Electrostatics, BuiltInPotentialZeroAtIntrinsic) {
    // Nd = Na = n_i -> V_bi = 0
    const double V_T = PhysicsEngine::thermalVoltage(300.0);
    EXPECT_NEAR(PhysicsEngine::builtInPotential(1.0e10, 1.0e10, 1.0e10, V_T),
                0.0, 1.0e-9);
}

TEST(Electrostatics, EquilibriumDensityAgreesAtZeroPsi) {
    // psi=0, phi=0 -> n = n_i
    const double V_T = PhysicsEngine::thermalVoltage(300.0);
    EXPECT_NEAR(
        PhysicsEngine::equilibriumElectronDensity(1.0e10, 0.0, V_T),
        1.0e10, 1.0e-6);
    EXPECT_NEAR(
        PhysicsEngine::equilibriumHoleDensity(1.0e10, 0.0, V_T),
        1.0e10, 1.0e-6);
}

TEST(Electrostatics, OhmicContactPsiAsymmetry) {
    // n-contact psi positive; p-contact psi negative (Sze Eq. 2.7).
    const double V_T = PhysicsEngine::thermalVoltage(300.0);
    const double n_i = 1.0e10;
    const double psi_n = PhysicsEngine::ohmicContactPsi(
        0.0, /*Nd*/1.0e17, /*Na*/1.0, n_i, V_T);
    const double psi_p = PhysicsEngine::ohmicContactPsi(
        0.0, /*Nd*/1.0, /*Na*/1.0e17, n_i, V_T);
    EXPECT_GT(psi_n,  0.0);
    EXPECT_LT(psi_p,  0.0);
    // Magnitudes match (same N): symmetric junction
    EXPECT_NEAR(psi_n, -psi_p, 1.0e-9);
}


// =============================================================================
// Capacitance helpers (Phase 4)
// =============================================================================
TEST(Capacitance, DepletionWidthShrinksUnderForwardBias) {
    using P = PhysicsEngine;
    const double V_bi = 0.7, eps_r = 11.7;
    const double W_eq = P::depletionWidthFlat(
        1.0e17, 1.0e17, V_bi, 0.0, eps_r);
    const double W_fwd = P::depletionWidthFlat(
        1.0e17, 1.0e17, V_bi, 0.4, eps_r);
    const double W_rev = P::depletionWidthFlat(
        1.0e17, 1.0e17, V_bi, -2.0, eps_r);
    EXPECT_LT(W_fwd, W_eq);
    EXPECT_GT(W_rev, W_eq);
}

TEST(Capacitance, DepletionCapacitanceScalesInverselyWithW) {
    using P = PhysicsEngine;
    // C = eps_s A / W -- doubling W halves C
    const double Cz = P::depletionCapacitanceFlat(
        1.0e17, 1.0e17, 0.7, 0.0, 11.7, 1.0);
    const double Cr = P::depletionCapacitanceFlat(
        1.0e17, 1.0e17, 0.7, -3.0, 11.7, 1.0);
    EXPECT_LT(Cr, Cz);
}

TEST(Capacitance, DiffusionCapacitanceProportionalToCurrent) {
    using P = PhysicsEngine;
    const double V_T = P::thermalVoltage(300.0);
    const double Cd1 = P::diffusionCapacitance(1.0e-3, 1.0e-7, V_T);
    const double Cd2 = P::diffusionCapacitance(2.0e-3, 1.0e-7, V_T);
    EXPECT_NEAR(Cd2 / Cd1, 2.0, 1.0e-6);
}


// =============================================================================
// Local mobility (Phase 4)  -- mu(N, E) recovers low-field limits
// =============================================================================
TEST(LocalMobility, ZeroFieldEqualsMatthiessenLowField) {
    const auto& mat = material::byKind(material::Kind::Silicon);
    const double T = 300.0, N = 1.0e16;
    const double mu_low = PhysicsEngine::matthiessenMobilityElectron(
        mat, T, N);
    const double mu_E0  = PhysicsEngine::localMobilityElectron(
        mat, T, N, 0.0);
    EXPECT_NEAR(mu_E0, mu_low, 1.0e-6 * mu_low);
}

TEST(LocalMobility, HighFieldSaturatesVelocity) {
    // v = mu_eff * E should plateau near v_sat at very high field.
    const auto& mat = material::byKind(material::Kind::Silicon);
    const double T = 300.0, N = 1.0e16;
    const double E = 5.0e5;   // 500 kV/cm -- deep in saturation
    const double mu_eff = PhysicsEngine::localMobilityElectron(mat, T, N, E);
    const double v      = mu_eff * E;
    EXPECT_LT(v, 1.3 * mat.v_sat_n);
}


// =============================================================================
// PN diode benchmark -- equilibrium V_bi from a painted junction
// -----------------------------------------------------------------------------
// Programmatic painter test: stamp Nd / Na on a DriftDiffusion grid, run
// equilibrium Poisson, and check that the psi step matches the textbook
// built-in potential within 5%.  Sze Sec. 2.2.
// =============================================================================
#include "DriftDiffusion.hpp"
TEST(PNJunction, EquilibriumVbiMatchesTextbook) {
    PhysicsEngine pe(material::Kind::Silicon);
    pe.setTemperature(300.0);
    pe.setDopingType(DopingType::Intrinsic);   // grid-defined doping
    DriftDiffusion dd(60, 16);
    dd.configureForMaterial(pe.getMaterial());
    dd.setCellPitchCm(1.0e-5f);    // 100 nm pitch -> 6 um device

    const double Nd = 1.0e17, Na = 1.0e17;
    for (int j = 0; j < dd.height(); ++j) {
        for (int i = 0; i < dd.width(); ++i) {
            if (i < dd.width() / 2) dd.setDopingAt(i, j, 0.0, Na);
            else                    dd.setDopingAt(i, j, Nd,  0.0);
        }
    }
    dd.setDeviceMode(DeviceMode::Painter);

    // Plenty of sweeps for tight equilibrium convergence.
    const double n_i = pe.getIntrinsicCarrier();
    const double V_T = PhysicsEngine::thermalVoltage(300.0);
    dd.solvePoisson(n_i, V_T, pe.getMaterial().epsilon_r, 5000, 0.85);

    // Measure psi step: right column (n-side) minus left column (p-side).
    const double psi_n = dd.psiAt(dd.width() - 2, dd.height() / 2);
    const double psi_p = dd.psiAt(1,              dd.height() / 2);
    const double psi_step = psi_n - psi_p;

    const double V_bi_textbook =
        PhysicsEngine::builtInPotential(Nd, Na, n_i, V_T);
    // Within 5% -- depletion-approximation tail at the contacts plus
    // discretisation; closer convergence comes from a finer mesh.
    EXPECT_NEAR(psi_step, V_bi_textbook, 0.05 * V_bi_textbook);
}
