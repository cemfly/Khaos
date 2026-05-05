// =============================================================================
// tests/test_physics.cpp
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
