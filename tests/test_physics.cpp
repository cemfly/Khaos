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
